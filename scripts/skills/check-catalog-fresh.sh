#!/bin/sh
# CI entrypoint (STRATEGY.md section 5): verify the generated skills tool
# catalog is fresh with respect to docs/mcp/tools/**. Exits non-zero on drift.
# Regenerate with: python3 scripts/skills/gen-catalog.py
exec python3 "$(dirname "$0")/gen-catalog.py" --check
