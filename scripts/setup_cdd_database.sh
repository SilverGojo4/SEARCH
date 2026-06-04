#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ===================== Usage =====================
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  echo
  echo "Example:"
  echo "  bash scripts/setup_cdd_database.sh ."
  exit 1
fi

# Convert project directory to absolute path
PROJECT_DIR=$(cd "$1" && pwd)

# ===================== Paths =====================
CDD_ARCHIVE="$PROJECT_DIR/data/external/cdd/archive/cdd.tar.gz"
CDD_PROFILE_DIR="$PROJECT_DIR/data/external/cdd/profiles"
CDD_RPSBLAST_DIR="$PROJECT_DIR/data/external/cdd/rpsblast"
CDD_METADATA_DIR="$PROJECT_DIR/data/external/cdd/metadata"
CDD_METADATA_TABLE="$CDD_METADATA_DIR/cddid.tbl"

CONDA_ENV="search"
DB_NAME="Cdd"
SMP_LIST="$CDD_RPSBLAST_DIR/${DB_NAME}.pn"

# ===================== Step 1: Check Files =====================
echo "🔍 Checking required files and directories ..."

if [ ! -f "$CDD_ARCHIVE" ]; then
  echo "❌ Missing CDD archive:"
  echo "  $CDD_ARCHIVE"
  echo
  echo "Please download cdd.tar.gz first and place it at:"
  echo "  data/external/cdd/archive/cdd.tar.gz"
  echo
  echo "Download source:"
  echo "  https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz"
  exit 1
fi

if [ ! -f "$CDD_METADATA_TABLE" ]; then
  echo "⚠️ CDD metadata table not found:"
  echo "  $CDD_METADATA_TABLE"
  echo
  echo "This is not required for makeprofiledb, but may be needed for downstream annotation."
  echo "You can prepare it with:"
  echo "  mkdir -p data/external/cdd/metadata"
  echo "  wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O data/external/cdd/metadata/cddid.tbl.gz"
  echo "  gunzip -k data/external/cdd/metadata/cddid.tbl.gz"
  echo
fi

mkdir -p "$CDD_PROFILE_DIR"
mkdir -p "$CDD_RPSBLAST_DIR"

# ===================== Step 2: Extract CDD Profiles =====================
echo "📦 Checking CDD profile extraction status ..."

if find "$CDD_PROFILE_DIR" -type f -name "*.smp" -print -quit | grep -q .; then
  echo "✅ Existing .smp files found. Skipping extraction."
  echo
else
  echo "📦 Extracting CDD profiles ..."
  echo "Input archive:"
  echo "  $CDD_ARCHIVE"
  echo "Output folder:"
  echo "  $CDD_PROFILE_DIR"
  echo

  # Use pv if available to show extraction progress.
  # Fall back to plain tar if pv is not installed.
  if command -v pv > /dev/null 2>&1; then
    pv "$CDD_ARCHIVE" | tar -xzf - -C "$CDD_PROFILE_DIR"
  else
    echo "⚠️ pv command not found. Extracting without progress bar."
    tar -xzf "$CDD_ARCHIVE" -C "$CDD_PROFILE_DIR"
  fi

  echo "✅ Extraction complete."
  echo
fi

# ===================== Step 3: Check .smp Files =====================
echo "🔍 Checking extracted .smp profile files ..."

# Do not use: ls "$CDD_PROFILE_DIR"/*.smp
# It may fail when there are too many .smp files.
find "$CDD_PROFILE_DIR" -type f -name "*.smp" | sort > "$SMP_LIST"

if [ ! -s "$SMP_LIST" ]; then
  echo "❌ No .smp files found in:"
  echo "  $CDD_PROFILE_DIR"
  echo
  echo "Please check whether cdd.tar.gz was extracted correctly."
  exit 1
fi

SMP_COUNT=$(wc -l < "$SMP_LIST")

echo "✅ Found $SMP_COUNT .smp files."
echo "📄 Profile list written to:"
echo "  $SMP_LIST"
echo

# ===================== Step 4: Activate Conda Environment =====================
echo "🐍 Activating Conda environment: $CONDA_ENV"

if ! command -v conda > /dev/null 2>&1; then
  echo "❌ conda command not found."
  echo "Please make sure Conda is installed and initialized."
  exit 1
fi

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# ===================== Step 5: Check makeprofiledb =====================
echo "🔍 Checking makeprofiledb availability ..."

if ! command -v makeprofiledb > /dev/null 2>&1; then
  echo "❌ makeprofiledb command not found in conda environment: $CONDA_ENV"
  echo
  echo "Please make sure BLAST+ is installed in this environment."
  echo "You can usually install it with:"
  echo "  conda install -c bioconda blast"
  exit 1
fi

# ===================== Step 6: Build RPS-BLAST Database =====================
echo "🔨 Building RPS-BLAST database ..."

cd "$CDD_RPSBLAST_DIR"

makeprofiledb \
  -in "$SMP_LIST" \
  -out "$DB_NAME" \
  -dbtype rps \
  -title "Conserved Domain Database (CDD)"

echo "✅ RPS-BLAST database build complete."
echo

# ===================== Step 7: Verify Outputs =====================
echo "📂 Verifying generated RPS-BLAST database files ..."

if ls "${DB_NAME}".* 1> /dev/null 2>&1; then
  ls -lh "${DB_NAME}".*
else
  echo "❌ No ${DB_NAME}.* files were generated."
  exit 1
fi

# ===================== Step 8: Deactivate Conda Environment =====================
conda deactivate

echo
echo "✅ CDD setup completed successfully!"
echo
echo "RPS-BLAST database path:"
echo "  $CDD_RPSBLAST_DIR/$DB_NAME"