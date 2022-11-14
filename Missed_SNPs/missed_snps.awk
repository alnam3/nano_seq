BEGIN {
    # Parameters.
    alternate_threshold = 10
    iupac_threshold = 0.3
    second_alternate_threshold = 1000

    # Codes.
    iupac["AC"] = "M"
    iupac["AG"] = "R"
    iupac["AT"] = "W"
    iupac["CG"] = "S"
    iupac["CT"] = "Y"
    iupac["GT"] = "K"
    # (+ symmetrical)
    for (pair in iupac) {
        split(pair, rev, "")
        iupac[rev[2] rev[1]] = iupac[pair]
    }

    seq_name = ""
}

# Output sequences fasta style.
seq_name != $1 {
    if (seq_name != "") {
        printf "\n"
    }
    seq_name=$1
    print ">", seq_name
}

# When no alternate is found, just spit out the ref base.
$4 == "." {
    printf $3
}

# Otherwise..
$4 != "." {

    # .. extract the counts and alternate bases found.
    ref = $3
    split($4, alternate, ",")
    split($5, counts, ",")

    # Rule out the situation with more than two credible alternates.
    if (counts[3] >= second_alternate_threshold) {
        print "Error: suspicious third count", \
              "in position", $2, "above." \
            | "cat 1>&2"
        exit 1
    }

    # Compare ref counts against major alternate.
    ratio = counts[2] / counts[1]

    if (ratio < iupac_threshold) {
        # Spit out the ref base again.
        printf ref
    } else if (ratio < alternate_threshold) {
        # The alternate is credible: pick the corresponding IUPAC code.
        printf iupac[ref alternate[1]]
    } else {
        # The alternate turns out to be the new reference.
        printf alternate[1]
    }
}

END {
    # That's it :)
    printf "\n"
}
