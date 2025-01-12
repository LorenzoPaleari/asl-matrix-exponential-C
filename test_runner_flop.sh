#!/bin/sh

# Set the range of n values
min_n=4
max_n=2048

# Initialize a flag for failed tests
tests_failed=0

# Ask the user for input
echo "Please enter a version name:"
read input

if echo "$input" | grep -Eq '^([[:alnum:]_]+)[[:space:]]+([[:alnum:]_]+)$'; then
    version=$(echo "$input" | awk '{print $1}')
    blas=$(echo "$input" | awk '{print $2}')
    echo "Version: $version"
    echo "n: $blas"
else
    echo "Input does not match the expected format."
fi

if [ "$blas" = "y" ]; then
  blas="S"
else
  blas="N"
fi

version_name=$version

# Run the 'make N={n}' command for each n value
rm -f ./data_collection/$version_name/flops.txt
while [ $min_n -le $max_n ]; do
  echo -n "Running test for N=$min_n "
  make cleanobj > /dev/null
  make N=$min_n V=$version_name FLOP=S ROOFLIN=N NO_BLA=$blas >> ./data_collection/$version_name/flops.txt
  exit_status=$?

  if [ $exit_status -ne 0 ]; then
    echo "❌"
    tests_failed=1
  else
    echo "✅"
  fi
  min_n=$((min_n * 2))
done

if [ $tests_failed -eq 1 ]; then
  exit 1
else
  exit 0
fi