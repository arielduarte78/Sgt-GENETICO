"""Pruebas de verificación de arquitectura e invariantes de dependencias vía AST.

Reglas del proyecto:
  - domain -> nada (solo biblioteca estándar y tipos puros; sin cvxpy, sin aleatoriedad).
  - evaluators -> domain. PROHIBIDO depender de search.
  - search -> evaluators y domain.
  - cases -> domain. PROHIBIDO depender de search o evaluators.
  - Dependencias en sentido contrario ESTRICTAMENTE PROHIBIDAS.
"""

from __future__ import annotations

import ast
from pathlib import Path


def _get_imports_from_file(file_path: Path) -> list[str]:
    """Extrae todos los módulos importados en runtime en un archivo .py mediante AST (ignorando TYPE_CHECKING)."""
    tree = ast.parse(file_path.read_text(encoding="utf-8"), filename=str(file_path))

    class RuntimeImportVisitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.imports: list[str] = []

        def visit_If(self, node: ast.If) -> None:
            # Ignorar bloques 'if TYPE_CHECKING:' (son tipos estáticos puros, sin dependencia en runtime)
            if isinstance(node.test, ast.Name) and node.test.id == "TYPE_CHECKING":
                return
            if isinstance(node.test, ast.Attribute) and node.test.attr == "TYPE_CHECKING":
                return
            self.generic_visit(node)

        def visit_Import(self, node: ast.Import) -> None:
            for alias in node.names:
                self.imports.append(alias.name)

        def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
            if node.module:
                self.imports.append(node.module)

    visitor = RuntimeImportVisitor()
    visitor.visit(tree)
    return visitor.imports



def test_domain_has_zero_external_solver_or_higher_layer_dependencies() -> None:
    """Verifica que el dominio matemático sea puro e independiente."""
    domain_dir = Path("src/kuramoto_templates/domain")
    py_files = list(domain_dir.glob("*.py"))
    assert py_files, "No se encontraron archivos en domain/"

    prohibited_prefixes = (
        "kuramoto_templates.evaluators",
        "kuramoto_templates.search",
        "kuramoto_templates.cases",
        "kuramoto_templates.io",
        "cvxpy",
        "clarabel",
        "scipy.optimize",
        "random",
    )

    for py_file in py_files:
        imports = _get_imports_from_file(py_file)
        for imp in imports:
            for prohibited in prohibited_prefixes:
                assert not imp.startswith(prohibited), (
                    f"Violación arquitectónica en {py_file.name}: "
                    f"el dominio puro no puede importar '{imp}' ({prohibited})."
                )


def test_evaluators_layer_does_not_depend_on_search() -> None:
    """Verifica que evaluators no importe el motor de búsqueda."""
    evaluators_dir = Path("src/kuramoto_templates/evaluators")
    py_files = list(evaluators_dir.glob("*.py"))
    assert py_files, "No se encontraron archivos en evaluators/"

    for py_file in py_files:
        imports = _get_imports_from_file(py_file)
        for imp in imports:
            assert not imp.startswith("kuramoto_templates.search"), (
                f"Violación arquitectónica en {py_file.name}: "
                f"evaluators no puede depender de search ('{imp}')."
            )


def test_cases_layer_does_not_depend_on_search_or_evaluators() -> None:
    """Verifica que cases no dependa de search ni de evaluators."""
    cases_dir = Path("src/kuramoto_templates/cases")
    py_files = list(cases_dir.glob("*.py"))
    assert py_files, "No se encontraron archivos en cases/"

    prohibited_prefixes = (
        "kuramoto_templates.search",
        "kuramoto_templates.evaluators",
    )

    for py_file in py_files:
        imports = _get_imports_from_file(py_file)
        for imp in imports:
            for prohibited in prohibited_prefixes:
                assert not imp.startswith(prohibited), (
                    f"Violación arquitectónica en {py_file.name}: "
                    f"cases no puede depender de '{prohibited}' (importa '{imp}')."
                )
