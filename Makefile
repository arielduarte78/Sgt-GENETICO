.PHONY: env test lint reproduce clean help

export UV_LINK_MODE := copy
UV ?= uv

help:
	@echo "Comandos disponibles:"
	@echo "  make env        - Crea el entorno virtual e instala las dependencias pineadas con uv"
	@echo "  make test       - Ejecuta la suite de pruebas unitarias y de arquitectura con pytest"
	@echo "  make lint       - Ejecuta ruff y mypy en modo estricto"
	@echo "  make reproduce  - Ejecuta las verificaciones de reproducibilidad del pipeline"
	@echo "  make clean      - Limpia caches y archivos temporales de Python"

env:
	$(UV) venv --python 3.12
	$(UV) pip install --link-mode=copy -e .[dev]

test:
	$(UV) run pytest tests

lint:
	$(UV) run ruff check src tests
	$(UV) run mypy --strict src tests

reproduce:
	@echo "[REPRODUCE] Verificando entorno e invariantes del paquete..."
	$(UV) run pytest tests
	@echo "[REPRODUCE] Estado del esqueleto: conforme a las especificaciones cientificas."

clean:
	rm -rf .pytest_cache .mypy_cache .ruff_cache build dist *.egg-info
