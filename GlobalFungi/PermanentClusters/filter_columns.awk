#!/usr/bin/awk -f

BEGIN {
    FS = "\t";  # Set input field separator to tab
    OFS = "\t"; # Set output field separator to tab

    # Set default threshold to 1000 if not provided
    threshold = 1000;
    if (ARGC > 1) {
        threshold = ARGV[1];
        ARGC--;  # Decrease ARGC to ignore the threshold parameter in input file processing
        ARGV[1] = ARGV[2];  # Shift ARGV to ignore the first argument (threshold) and use the input file
    }
}

# This block runs on the header line to store header names
NR == 1 {
    for (i = 2; i <= NF; i++) {
        header[i] = $i;
    }
    next;
}

{
    # Calculate sum for each column and store the rows
    for (i = 2; i <= NF; i++) {
        column_sum[i] += $i;
    }

    # Store each row for later output
    data[NR] = $0;
}

# END block to print selected columns and filter rows
END {
    # Select columns based on the sum compared to the threshold
    selected_columns[1] = 1;  # Always keep the first column (cluster_name)
    for (i = 2; i <= NF; i++) {
        if (column_sum[i] >= threshold) {
            selected_columns[i] = 1;
        }
    }

    # Print header for selected columns
    printf "%s", "cluster_name";
    for (i = 2; i <= NF; i++) {
        if (selected_columns[i]) {
            printf "%s%s", OFS, header[i];
        }
    }
    printf "\n";

    # Process each row again to print rows where sum > 0
    for (row = 2; row <= NR; row++) {
        split(data[row], fields, FS);
        row_sum = 0;

        # Calculate the row sum
        for (i = 2; i <= NF; i++) {
            row_sum += fields[i];
        }

        # Only print the row if sum > 0
        if (row_sum > 0) {
            printf "%s", fields[1];  # Print first field (cluster_name)
            for (i = 2; i <= NF; i++) {
                if (selected_columns[i]) {
                    printf "%s%s", OFS, fields[i];
                }
            }
            printf "\n";
        }
    }
}
