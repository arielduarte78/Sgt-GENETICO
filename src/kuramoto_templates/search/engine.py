"""Island-model Genetic Algorithm search engine with ring migration and deterministic checkpoints."""

from __future__ import annotations

from collections.abc import Callable
import json
from pathlib import Path
import time
from typing import Any
import numpy as np

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.results import EvalResult
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.search.checkpoint import (
    CheckpointState,
    load_checkpoint,
    save_checkpoint_atomic,
)
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.genome import (
    Genome,
    canonicalize,
    to_phenotype,
)
from kuramoto_templates.search.manifest import create_run_manifest, save_manifest
from kuramoto_templates.search.operators.dispatcher import (
    OperatorRegistry,
    select_and_apply_operator,
)
from kuramoto_templates.search.runlog import RunLogger
from kuramoto_templates.search.seeds import (
    create_c4_seed,
    create_s7_direct_seed,
)
from kuramoto_templates.search.selection import (
    Individual,
    compare_individuals,
    compute_f_sel,
    deterministic_crowding_replace,
    tournament_selection,
)

EvaluatorFunc = Callable[[tuple[float, ...], tuple[float, ...]], EvalResult]


def _individual_to_dict(ind: Individual) -> dict[str, Any]:
    """Serialize Individual to a JSON-compatible dictionary."""
    return {
        "id": ind.id,
        "genome": {
            "m": ind.genome.m,
            "K_ref": ind.genome.K_ref,
            "k": list(ind.genome.k),
            "sgn": list(ind.genome.sgn),
            "ell": list(ind.genome.ell),
            "u": list(ind.genome.u),
        },
        "phenotype_x": list(ind.phenotype_x),
        "phenotype_phi": list(ind.phenotype_phi),
        "s_star": ind.s_star,
        "lam_star": ind.lam_star,
        "lam_normalized": ind.lam_normalized,
        "f": ind.f,
        "f_sel": ind.f_sel,
        "parent_ids": list(ind.parent_ids),
        "operator_name": ind.operator_name,
        "generation": ind.generation,
        "solver_status": ind.solver_status,
        "eval_time": ind.eval_time,
    }


def _dict_to_individual(d: dict[str, Any]) -> Individual:
    """Deserialize Individual from a JSON-compatible dictionary."""
    g_dict = d["genome"]
    genome = Genome(
        m=int(g_dict["m"]),
        K_ref=int(g_dict["K_ref"]),
        k=tuple(int(v) for v in g_dict["k"]),
        sgn=tuple(int(v) for v in g_dict["sgn"]),
        ell=tuple(float(v) for v in g_dict["ell"]),
        u=tuple(float(v) for v in g_dict["u"]),
    )
    return Individual(
        id=int(d["id"]),
        genome=genome,
        phenotype_x=tuple(float(v) for v in d["phenotype_x"]),
        phenotype_phi=tuple(float(v) for v in d["phenotype_phi"]),
        s_star=float(d["s_star"]),
        lam_star=float(d["lam_star"]),
        lam_normalized=float(d["lam_normalized"]),
        f=float(d["f"]),
        f_sel=float(d["f_sel"]),
        parent_ids=tuple(int(v) for v in d["parent_ids"]),
        operator_name=str(d["operator_name"]),
        generation=int(d.get("generation", 0)),
        solver_status=str(d.get("solver_status", "OPTIMAL")),
        eval_time=float(d.get("eval_time", 0.0)),
    )


class GAEngine:
    """Genetic Algorithm Engine orchestrating island populations, SDP evaluations, and logging."""

    def __init__(
        self,
        config: GAConfig,
        output_dir: Path | str,
        evaluator_fn: EvaluatorFunc | None = None,
        run_id: str | None = None,
        allow_dirty: bool = False,
        no_seeds: bool = False,
        no_split_merge: bool = False,
    ) -> None:
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.run_id = run_id or f"ga_{int(time.time())}"
        self.allow_dirty = allow_dirty
        self.no_seeds = no_seeds or not config.use_seeds
        self.no_split_merge = no_split_merge or not config.allow_split_merge

        self.evaluator_fn = evaluator_fn
        if self.evaluator_fn is None:
            sdp_eval = SDPEvaluator(EvalConfig(mu_t=config.mu_t, solver_tol=config.tol))
            self._eval_func: EvaluatorFunc = sdp_eval.evaluate_phase1
        else:
            self._eval_func = self.evaluator_fn

        self.logger = RunLogger(self.output_dir / "run.jsonl", run_id=self.run_id)
        self.checkpoint_dir = self.output_dir / "checkpoints"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

        self.cache: dict[str, dict[str, Any]] = {}
        self.cache_hits: int = 0
        self.next_id: int = 1

        self.islands: list[list[Individual]] = [[] for _ in range(config.num_islands)]
        self.hall_of_fame_global: list[Individual] = []
        self.hall_of_fame_per_m: dict[int, Individual] = {}

        self.registry = OperatorRegistry(
            probabilities=dict(config.operator_probs),
            allow_split_merge=not self.no_split_merge,
        )

        # SeedSequence initialized once deterministically with invariant pool size
        self.seed_seq = np.random.SeedSequence(config.seed)
        self.spawned_seeds = [
            int(s.generate_state(1)[0])
            for s in self.seed_seq.spawn(10000)
        ]

        # Save run manifest
        manifest = create_run_manifest(
            config=config,
            run_id=self.run_id,
            allow_dirty=self.allow_dirty,
            no_seeds=self.no_seeds,
            no_split_merge=self.no_split_merge,
            spawned_seeds=self.spawned_seeds,
        )
        save_manifest(manifest, self.output_dir / "manifest.json")

        self.current_generation: int = 0

    def _get_island_seed(self, gen: int, island_idx: int) -> int:
        """Deterministic seed lookup for (generation, island)."""
        idx = gen * self.config.num_islands + island_idx
        if idx < len(self.spawned_seeds):
            return self.spawned_seeds[idx]
        return int(self.config.seed + idx * 7919)

    def _evaluate(
        self,
        genome: Genome,
        parent_ids: tuple[int, ...],
        operator_name: str,
        island: int,
        generation: int,
        phenotype_override: tuple[tuple[float, ...], tuple[float, ...]] | None = None,
    ) -> Individual:
        """Evaluate a genome or phenotype override with caching and JSONL logging."""
        if phenotype_override is not None:
            x, phi = phenotype_override
        else:
            x, phi = to_phenotype(genome)

        cache_key = json.dumps(
            [round(float(v), 10) for v in x] + [round(float(v), 10) for v in phi]
        )

        if cache_key in self.cache:
            self.cache_hits += 1
            entry = self.cache[cache_key]
            s_star = float(entry["s_star"])
            lam_star = float(entry["lam_star"])
            lam_normalized = float(entry["lam_normalized"])
            solver_status = str(entry["solver_status"])
            eval_time = float(entry["eval_time"])
        else:
            t0 = time.perf_counter()
            res = self._eval_func(x, phi)
            eval_time = time.perf_counter() - t0

            s_star = float(res.s) if res.s is not None else 0.0
            lam_star = float(res.lam) if res.lam is not None else 0.0
            lam_normalized = (
                float(res.lam_normalized) if res.lam_normalized is not None else 0.0
            )
            solver_status = str(res.status)

            self.cache[cache_key] = {
                "s_star": s_star,
                "lam_star": lam_star,
                "lam_normalized": lam_normalized,
                "solver_status": solver_status,
                "eval_time": eval_time,
            }

        f = compute_continuous_fitness(lam=lam_star, s=s_star, tol=self.config.tol)
        f_sel = compute_f_sel(s_star, lam_normalized, tol=self.config.tol)

        ind = Individual(
            id=self.next_id,
            genome=genome,
            phenotype_x=x,
            phenotype_phi=phi,
            s_star=s_star,
            lam_star=lam_star,
            lam_normalized=lam_normalized,
            f=f,
            f_sel=f_sel,
            parent_ids=parent_ids,
            operator_name=operator_name,
            generation=generation,
            solver_status=solver_status,
            eval_time=eval_time,
        )
        self.next_id += 1

        self.logger.log_evaluation(
            eval_id=ind.id,
            parent_ids=parent_ids,
            operator=operator_name,
            f=f,
            s_star=s_star,
            lam_raw=lam_star,
            lam_norm=lam_normalized,
            solver_status=solver_status,
            eval_time=eval_time,
            m=genome.m,
            island=island,
            generation=generation,
        )

        self._update_hall_of_fame(ind)
        return ind

    def _update_hall_of_fame(self, ind: Individual) -> None:
        """Update global and per-m Hall of Fame archives."""
        # Global Hall of Fame
        if not self.hall_of_fame_global:
            self.hall_of_fame_global.append(ind)
        else:
            best = self.hall_of_fame_global[0]
            if compare_individuals(ind, best, self.config.fitness_noise) > 0:
                self.hall_of_fame_global.insert(0, ind)
                if len(self.hall_of_fame_global) > 20:
                    self.hall_of_fame_global.pop()

        # Per-m Hall of Fame
        m = ind.genome.m
        if m not in self.hall_of_fame_per_m:
            self.hall_of_fame_per_m[m] = ind
        else:
            curr_best = self.hall_of_fame_per_m[m]
            if compare_individuals(ind, curr_best, self.config.fitness_noise) > 0:
                self.hall_of_fame_per_m[m] = ind

    def initialize_population(self) -> None:
        """Create and evaluate the initial population across all islands."""
        init_rng = np.random.default_rng(self.config.seed)

        for island_idx in range(self.config.num_islands):
            pop: list[Individual] = []

            # Inject seeds if enabled
            if not self.no_seeds:
                s7_seed = create_s7_direct_seed()
                ind_s7 = self._evaluate(
                    genome=s7_seed.genome,
                    parent_ids=(),
                    operator_name="seed_s7",
                    island=island_idx,
                    generation=0,
                    phenotype_override=(s7_seed.phenotype_x, s7_seed.phenotype_phi),
                )
                pop.append(ind_s7)

                c4_seed = create_c4_seed()
                ind_c4 = self._evaluate(
                    genome=c4_seed.genome,
                    parent_ids=(),
                    operator_name="seed_c4",
                    island=island_idx,
                    generation=0,
                )
                pop.append(ind_c4)

            # Fill the rest with randomized valid genomes
            while len(pop) < self.config.island_pop_size:
                m = int(init_rng.integers(2, 11))
                kref = int(init_rng.choice([4, 8, 12, 16]))
                k = tuple(int(init_rng.integers(0, kref)) for _ in range(m))
                sgn = tuple(int(init_rng.choice([-1, 1])) for _ in range(m))
                ell = tuple(float(init_rng.uniform(-6.0, -0.5)) for _ in range(m))
                u = tuple(float(init_rng.normal(0.0, 1.0)) for _ in range(m))

                raw_genome = Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u)
                canonical_genome = canonicalize(raw_genome)

                ind = self._evaluate(
                    genome=canonical_genome,
                    parent_ids=(),
                    operator_name="random_init",
                    island=island_idx,
                    generation=0,
                )
                pop.append(ind)

            self.islands[island_idx] = pop

    def _ring_migration(self) -> None:
        """Migrate best individuals in a directional ring between islands."""
        if self.config.num_islands <= 1:
            return

        migrants: list[list[Individual]] = []
        for pop in self.islands:
            # Sort island descending by fitness
            sorted_pop = sorted(
                pop,
                key=lambda ind: ind.f_sel,
                reverse=True,
            )
            migrants.append(sorted_pop[: self.config.migration_count])

        for i in range(self.config.num_islands):
            target_idx = (i + 1) % self.config.num_islands
            immigrants = migrants[i]
            target_pop = self.islands[target_idx]

            # Replace weakest individuals in target island if immigrants are better
            target_pop.sort(key=lambda ind: ind.f_sel)
            for j, imm in enumerate(immigrants):
                if (
                    compare_individuals(
                        imm, target_pop[j], self.config.fitness_noise
                    )
                    > 0
                ):
                    target_pop[j] = imm

    def run(self) -> dict[str, Any]:
        """Execute evolutionary search across generations."""
        if self.current_generation == 0 and not any(self.islands):
            self.initialize_population()

        start_gen = self.current_generation + 1
        total_gens = self.config.num_generations

        for gen in range(start_gen, total_gens + 1):
            self.current_generation = gen

            for island_idx in range(self.config.num_islands):
                island_rng = np.random.default_rng(
                    self._get_island_seed(gen, island_idx)
                )
                current_pop = self.islands[island_idx]
                next_pop: list[Individual] = []

                # Keep elite individual (mu+lambda preservation)
                elite = max(
                    current_pop,
                    key=lambda ind: ind.f_sel,
                )
                next_pop.append(elite)

                while len(next_pop) < self.config.island_pop_size:
                    # Tournament selection for parents
                    p1 = tournament_selection(
                        current_pop,
                        k=self.config.tournament_k,
                        fitness_noise=self.config.fitness_noise,
                        rng=island_rng,
                    )
                    p2 = tournament_selection(
                        current_pop,
                        k=self.config.tournament_k,
                        fitness_noise=self.config.fitness_noise,
                        rng=island_rng,
                    )

                    # Apply operator
                    child_genome, op_name = select_and_apply_operator(
                        p1.genome,
                        p2.genome,
                        self.registry,
                        island_rng,
                    )
                    child_canonical = canonicalize(child_genome)

                    child_ind = self._evaluate(
                        genome=child_canonical,
                        parent_ids=(p1.id, p2.id),
                        operator_name=op_name,
                        island=island_idx,
                        generation=gen,
                    )

                    # Deterministic crowding replacement against nearest parent
                    next_winner, _ = deterministic_crowding_replace(
                        p1, p2, child_ind, p2, self.config.fitness_noise
                    )
                    next_pop.append(next_winner)

                self.islands[island_idx] = next_pop

                # Log generation metrics
                best_ind = max(next_pop, key=lambda ind: ind.f_sel)
                mean_f = float(np.mean([ind.f for ind in next_pop]))
                self.logger.log_generation_summary(
                    island=island_idx,
                    generation=gen,
                    best_f=best_ind.f,
                    best_m=best_ind.genome.m,
                    best_lam_raw=best_ind.lam_star,
                    mean_f=mean_f,
                    cache_hits=self.cache_hits,
                    evals_count=len(self.cache),
                )

            # Ring migration
            if gen % self.config.migration_interval == 0:
                self._ring_migration()

            # Atomic checkpointing
            if (
                gen % self.config.checkpoint_interval == 0
                or gen == total_gens
            ):
                self.save_checkpoint(gen)

        best_global = (
            self.hall_of_fame_global[0] if self.hall_of_fame_global else None
        )
        return {
            "run_id": self.run_id,
            "generations": total_gens,
            "best_f": best_global.f if best_global else -1.0,
            "best_m": best_global.genome.m if best_global else -1,
            "best_lam": best_global.lam_star if best_global else 0.0,
            "best_s": best_global.s_star if best_global else 1.0,
            "cache_hits": self.cache_hits,
            "total_evaluations": len(self.cache),
            "hall_of_fame_per_m": {
                m: _individual_to_dict(ind)
                for m, ind in self.hall_of_fame_per_m.items()
            },
        }

    def save_checkpoint(self, gen: int) -> Path:
        """Create an atomic checkpoint at generation gen."""
        state = CheckpointState(
            run_id=self.run_id,
            generation=gen,
            island_populations=[
                [_individual_to_dict(ind) for ind in pop] for pop in self.islands
            ],
            hall_of_fame_global=[
                _individual_to_dict(ind) for ind in self.hall_of_fame_global
            ],
            hall_of_fame_per_m={
                m: _individual_to_dict(ind)
                for m, ind in self.hall_of_fame_per_m.items()
            },
            cache=self.cache,
            cache_hits=self.cache_hits,
            next_individual_id=self.next_id,
        )
        cp_path = self.checkpoint_dir / f"checkpoint_gen_{gen}.json"
        save_checkpoint_atomic(state, cp_path)
        return cp_path

    @classmethod
    def resume(
        cls,
        checkpoint_path: Path | str,
        total_generations: int | None = None,
        evaluator_fn: EvaluatorFunc | None = None,
    ) -> GAEngine:
        """Resume search run exactly from an existing checkpoint file."""
        state = load_checkpoint(checkpoint_path)
        cp_path = Path(checkpoint_path)
        output_dir = cp_path.parent.parent

        manifest_path = output_dir / "manifest.json"
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)

        raw_config = manifest["config"]
        if total_generations is not None:
            raw_config["num_generations"] = total_generations

        config = GAConfig(**raw_config)

        engine = cls(
            config=config,
            output_dir=output_dir,
            evaluator_fn=evaluator_fn,
            run_id=state.run_id,
            allow_dirty=manifest.get("is_dirty", True),
            no_seeds=manifest["ablation_flags"].get("no_seeds", False),
            no_split_merge=manifest["ablation_flags"].get("no_split_merge", False),
        )

        engine.current_generation = state.generation
        engine.cache = dict(state.cache)
        engine.cache_hits = state.cache_hits
        engine.next_id = state.next_individual_id

        engine.islands = [
            [_dict_to_individual(ind_d) for ind_d in pop]
            for pop in state.island_populations
        ]
        engine.hall_of_fame_global = [
            _dict_to_individual(ind_d) for ind_d in state.hall_of_fame_global
        ]
        engine.hall_of_fame_per_m = {
            m: _dict_to_individual(ind_d)
            for m, ind_d in state.hall_of_fame_per_m.items()
        }

        return engine
