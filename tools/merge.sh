# original implementation code base without tests
find . \
    \( -path './Molmodel/examplesNotReadyYet' \
    -o -path './Molmodel/examples' \
    -o -path './Molmodel/rna-dynamics' \
    -o -path './pybind11' \
    -o -path './openmm' \
    -o -path './tools' \
    -o -path './tests' \) \
    -prune -o \
    -type f \
    \( -name '*.py' \
    -o -name '*.c' \
    -o -name '*.cpp' \
    -o -name '*.h' \
    -o -name '*.hpp' \) \
    -print |
sort |
while IFS= read -r f; do
    echo "<file path=\"$f\">"
    cat "$f"
    echo
    echo "</file>"
    echo
done > codebase_original_with_tests.xml


(find python/robosample -type f -name '*.py'; find src -type f -name '*.cpp'; find include -type f -name '*.hpp'; find tests -type f -name '*.hpp'; find tests -type f -name '*.cpp') | sort | while IFS= read -r f; do echo "<file path=\"$f\">"; cat "$f"; echo; echo "</file>"; echo; done > llm/codebase_current.xml