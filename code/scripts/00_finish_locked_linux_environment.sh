#!/usr/bin/env bash
# Repair prefix-local GUI caches after an explicit Conda lock installation.
# Bioconductor data post-link ordering is handled by the R source installer.
set -euo pipefail
[[ $# == 1 ]] || { printf 'Usage: %s /absolute/conda/prefix\n' "$0" >&2; exit 2; }
AIS_DEPENDENCY_PREFIX="$(cd "$1" && pwd)"
[[ "$(uname -s)" == Linux ]] || { printf 'This helper is for Linux.\n' >&2; exit 2; }
export PATH="$AIS_DEPENDENCY_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$AIS_DEPENDENCY_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
if [[ -x "$AIS_DEPENDENCY_PREFIX/bin/gdk-pixbuf-query-loaders" ]]; then
  "$AIS_DEPENDENCY_PREFIX/bin/gdk-pixbuf-query-loaders" > "$AIS_DEPENDENCY_PREFIX/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache"
fi
if [[ -x "$AIS_DEPENDENCY_PREFIX/bin/glib-compile-schemas" ]]; then
  "$AIS_DEPENDENCY_PREFIX/bin/glib-compile-schemas" "$AIS_DEPENDENCY_PREFIX/share/glib-2.0/schemas"
  "$AIS_DEPENDENCY_PREFIX/bin/gtk-update-icon-cache" -f -t "$AIS_DEPENDENCY_PREFIX/share/icons/hicolor"
  "$AIS_DEPENDENCY_PREFIX/bin/gtk-query-immodules-3.0" --update-cache
fi
printf 'Prefix-local Linux cache setup PASS. Run the locked R source installer next.\n'
