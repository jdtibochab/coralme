#!/usr/bin/env bash

# Embed the activation function
activate-python3-venv() {
    source /opt/python3-venv/$1/bin/activate
}

# Virtual environments
VENVS=(
  "p39-coralME"
  "p310-coralME"
  "p311-coralME"
  "p312-coralME"
  "p313-coralME"
  "p314-coralME"
)

# Display names for the columns
PY_NAMES=(
  "Python 3.9"
  "Python 3.10"
  "Python 3.11"
  "Python 3.12"
  "Python 3.13"
  "Python 3.14"
)

# Packages
PACKAGES=(
  "cobra"
  "numpy"
  "scipy"
  "pandas"
  "biopython"
  "sympy"
  "pint"
  "anyconfig"
  "requests"
)

# Gather versions
declare -A VERSIONS

for i in "${!VENVS[@]}"; do
    VENV="${VENVS[$i]}"
    if ! activate-python3-venv "$VENV" 2>/dev/null; then
        for PKG in "${PACKAGES[@]}"; do
            VERSIONS["$PKG,$i"]="missing"
        done
        continue
    fi

    FREEZE=$(pip freeze)

    for PKG in "${PACKAGES[@]}"; do
        VERSION=$(echo "$FREEZE" | grep -i "^${PKG}==" | cut -d= -f3)
        if [ -z "$VERSION" ]; then
            VERSION="not installed"
        fi
        VERSIONS["$PKG,$i"]="$VERSION"
    done
done

# Compute column widths
COL_WIDTHS=()
# Compute width of first column
max_len=0
for PKG in "${PACKAGES[@]}"; do
    (( ${#PKG} > max_len )) && max_len=${#PKG}
done
header="Package"
(( ${#header} > max_len )) && max_len=${#header}
COL_WIDTHS+=( $max_len )

# Other columns
for i in "${!PY_NAMES[@]}"; do
    max_len=${#PY_NAMES[$i]}
    for PKG in "${PACKAGES[@]}"; do
        V=${VERSIONS["$PKG,$i"]}
        (( ${#V} > max_len )) && max_len=${#V}
    done
    COL_WIDTHS+=( $max_len )
done

# Print header
printf "| %-*s " "${COL_WIDTHS[0]}" "Package"
for i in "${!PY_NAMES[@]}"; do
    printf "| %-*s " "${COL_WIDTHS[$((i+1))]}" "${PY_NAMES[$i]}"
done
echo "|"

# Print separator
printf "|"
for W in "${COL_WIDTHS[@]}"; do
    printf '%*s' "$W" '' | tr ' ' '-'
    printf "|"
done
echo

# Print each row
for PKG in "${PACKAGES[@]}"; do
    printf "| %-*s " "${COL_WIDTHS[0]}" "$PKG"
    for i in "${!VENVS[@]}"; do
        V=${VERSIONS["$PKG,$i"]}
        printf "| %-*s " "${COL_WIDTHS[$((i+1))]}" "$V"
    done
    echo "|"
done
