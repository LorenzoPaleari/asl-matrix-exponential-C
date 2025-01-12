#!/bin/sh

# Set the range of n values
min_n=2
max_n=16

# Initialize a flag for failed tests
tests_failed=0

# Ask the user for input
echo "Please enter a version name:"
read version_name

# Run the 'make N={n}' command for each n value
make cleanobj > /dev/null
for n in $(seq $min_n $max_n); do
  echo -n "Running test for N=$n "
  make N=$n V=base FLOP=N CALLGRIN=N >/dev/null 2>&1
  exit_status=$?

  if [ $exit_status -ne 0 ]; then
    echo "❌"
    tests_failed=1
  else
    echo "✅"
  fi
done

if [ $tests_failed -eq 1 ]; then
  exit 1
else
  exit 0
fi