base="/home/victor/Downloads/Robosample-5a2d451d38ca3c4a48f2d8c55f1c731e40459d64/python/robosample"
output="robosample_python.xml"

: > "$output"

find "$base" -type f -name "*.py" | sort | while IFS= read -r f; do
    echo "<file path=\"$f\">" >> "$output"
    cat "$f" >> "$output"
    echo -e "\n</file>\n" >> "$output"
done