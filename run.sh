find Inputs_and_Results/SKETCHES/anandhu -type f -name "*points.obj" | while read file; do
    echo "Processing: $file"
    # your_command "$file"
    ./build/example $file 128 10
done