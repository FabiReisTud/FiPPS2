#!/bin/bash
# combine_sources.sh
# Combines all files from include/ and source/ into a time-stamped combined file:
#   combined_2026-09-30_09-23-00_UTC.txt
# with a "===== filename =====" header before each file.

# --- Timecode (sortable, filename-safe: no spaces/colons) -------------------
TIMECODE="$(date -u '+%Y-%m-%d_%H-%M-%S_UTC')"
OUT="combined_${TIMECODE}.txt"
# ----------------------------------------------------------------------------

# Start with an empty output file
> "$OUT"

# Optional: keep a small header inside the file too (harmless, helps grepping)
echo "===== COMBINED_FILE_INFO =====" >> "$OUT"
echo "timecode: $TIMECODE"            >> "$OUT"
echo "generated_by: combine_sources.sh" >> "$OUT"
echo ""                              >> "$OUT"

# Loop over both folders
for dir in include source; do
    # Skip folder if it doesn't exist
    [ -d "$dir" ] || { echo "Warning: folder '$dir' not found, skipping." >&2; continue; }

    # find guarantees all files are found, even in subfolders;
    # sort for a deterministic order
    while IFS= read -r -d '' f; do
        # Skip the output file itself (in case it lives in one of these dirs)
        [ "$(realpath "$f")" = "$(realpath "$OUT")" ] && continue

        echo "===== $f =====" >> "$OUT"
        cat "$f" >> "$OUT"
        echo "" >> "$OUT"   # blank line so files never run into each other
    done < <(find "$dir" -type f -print0 | sort -z)
done

echo "Done. Wrote $(wc -l < "$OUT") lines to $OUT"
