"""Unit tests for GAConfig loading, validation, and SHA-256 calculation."""

import tempfile
from pathlib import Path

import pytest

from kuramoto_templates.search.config import GAConfig


def test_gaconfig_defaults() -> None:
    """GAConfig should have sensible defaults with required mu_t."""
    cfg = GAConfig(mu_t=0.690)
    assert cfg.mu_t == 0.690
    assert cfg.fitness_noise == 1e-9
    assert cfg.tol == 1e-8
    assert cfg.num_islands >= 1
    assert cfg.island_pop_size >= 2
    assert cfg.num_generations >= 1
    assert cfg.use_seeds is True
    assert cfg.allow_split_merge is True


def test_gaconfig_yaml_roundtrip_and_sha256() -> None:
    """YAML serialization/deserialization and sha256 checksum."""
    cfg = GAConfig(
        mu_t=0.6929077,
        fitness_noise=1e-9,
        num_islands=2,
        island_pop_size=10,
        num_generations=5,
        use_seeds=False,
        allow_split_merge=False,
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        yaml_path = Path(tmp_dir) / "config.yaml"
        cfg.to_yaml(yaml_path)

        loaded_cfg = GAConfig.from_yaml(yaml_path)
        assert loaded_cfg.mu_t == pytest.approx(0.6929077)
        assert loaded_cfg.use_seeds is False
        assert loaded_cfg.allow_split_merge is False

        sha1 = cfg.compute_sha256()
        sha2 = loaded_cfg.compute_sha256()
        assert len(sha1) == 64
        assert sha1 == sha2


def test_gaconfig_rejects_negative_noise_or_tol() -> None:
    """Config validation rejects negative fitness noise or invalid tolerances."""
    with pytest.raises(ValueError, match="fitness_noise"):
        GAConfig(mu_t=0.690, fitness_noise=-1e-9)

    with pytest.raises(ValueError, match="tol"):
        GAConfig(mu_t=0.690, tol=-1e-8)

    with pytest.raises(ValueError, match="h_eff_min"):
        GAConfig(mu_t=0.690, h_eff_min=-1e-4)



def test_gaconfig_min_max_m_bounds() -> None:
    """GAConfig should allow specifying min_m and max_m bounds and validate them."""
    cfg = GAConfig(mu_t=0.690, min_m=5, max_m=8)
    assert cfg.min_m == 5
    assert cfg.max_m == 8

    with pytest.raises(ValueError, match="min_m"):
        GAConfig(mu_t=0.690, min_m=1, max_m=8)

    with pytest.raises(ValueError, match="max_m"):
        GAConfig(mu_t=0.690, min_m=8, max_m=5)

