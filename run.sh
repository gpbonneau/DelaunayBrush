find flowrep -type f -name "*normal*" | while read file; do
    echo "Processing: $file"
    # your_command "$file"
    /build/example
done