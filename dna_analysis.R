# Programmer: Nita Hegde
# Purpose: Read a DNA sequence, create its reverse complement,
# and count nucleotides within each kilobase.
# Uses base R only.

# ------------------------------------------------------------
# Part 1: Read the complete DNA sequence
# ------------------------------------------------------------

fasta_lines <- readLines("chr1_GL383518v1_alt.fa")

# This script expects a FASTA file containing one sequence.
if (sum(grepl("^>", fasta_lines)) != 1) {
  stop("The FASTA file must contain exactly one sequence.")
}

# Remove the FASTA header, which starts with >.
sequence_lines <- fasta_lines[!grepl("^>", fasta_lines)]

# Join all sequence lines into one string.
dna <- paste(sequence_lines, collapse = "")

# Remove whitespace and convert lowercase bases to uppercase.
dna <- toupper(gsub("[[:space:]]", "", dna))

if (nchar(dna) < 800) {
  stop("The sequence must contain at least 800 bases.")
}

cat("Sequence length:", nchar(dna), "\n")

# R uses positions starting at 1.
cat("Part 1: 10th letter:", substr(dna, 10, 10), "\n")
cat("Part 1: 758th letter:", substr(dna, 758, 758), "\n")


# ------------------------------------------------------------
# Part 2: Create the reverse complement
# ------------------------------------------------------------

# Split the DNA string into individual letters.
bases <- strsplit(dna, split = "")[[1]]

# Reverse their order.
reversed_bases <- rev(bases)

# Join the reversed letters into one string.
reversed_dna <- paste(reversed_bases, collapse = "")

# Complement A <-> T and C <-> G.
# Also handle standard ambiguous DNA symbols.
# N represents an unknown base and remains N.
reverse_complement <- chartr(
  "ACGTRYSWKMBDHVN",
  "TGCAYRSWMKVHDBN",
  reversed_dna
)

cat(
  "Part 2: 79th letter:",
  substr(reverse_complement, 79, 79),
  "\n"
)

cat(
  "Part 2: Letters 500 through 800:",
  substr(reverse_complement, 500, 800),
  "\n"
)


# ------------------------------------------------------------
# Part 3: Count every letter within each kilobase
# ------------------------------------------------------------

# Each complete kilobase contains 1000 sequence positions.
# ceiling() includes a final block shorter than 1000 bases.
number_of_kb <- ceiling(length(bases) / 1000)

# Create an empty outer list with one element per kilobase.
kb_counts <- vector("list", number_of_kb)

# Include every observed letter and always include A, C, G, T.
letters_to_count <- sort(unique(c("A", "C", "G", "T", bases)))

for (i in seq_len(number_of_kb)) {
  
  start_position <- (i - 1) * 1000 + 1
  end_position <- min(i * 1000, length(bases))
  
  kb_bases <- bases[start_position:end_position]
  
  # A named list acts like a dictionary:
  # each letter is associated with its count.
  counts <- vector("list", length(letters_to_count))
  names(counts) <- letters_to_count
  
  for (letter in letters_to_count) {
    counts[[letter]] <- sum(kb_bases == letter)
  }
  
  kb_counts[[i]] <- counts
}

names(kb_counts) <- paste0("kb_", seq_len(number_of_kb))

# Display the counts for the first kilobase.
cat("\nPart 3: First kilobase counts:\n")
print(kb_counts[[1]])


# ------------------------------------------------------------
# Part 4: Create a four-column data frame from the list
# ------------------------------------------------------------

# Each row represents one kilobase.
# The four columns contain counts for A, C, G, and T.
nucleotide_df <- data.frame(
  A = integer(number_of_kb),
  C = integer(number_of_kb),
  G = integer(number_of_kb),
  T = integer(number_of_kb)
)

rownames(nucleotide_df) <- names(kb_counts)

for (i in seq_len(number_of_kb)) {
  
  counts <- kb_counts[[i]]
  
  nucleotide_df[i, ] <- c(
    counts[["A"]],
    counts[["C"]],
    counts[["G"]],
    counts[["T"]]
  )
}

cat("\nPart 4: Nucleotide counts per kilobase:\n")
print(nucleotide_df)

# Calculate row sums separately to retain exactly four columns.
row_totals <- rowSums(nucleotide_df)

cat("\nRow sums:\n")
print(row_totals)

# Identify rows whose A + C + G + T count is not 1000.
different_rows <- which(row_totals != 1000)

if (length(different_rows) == 0) {
  cat("\nAll rows have a sum of 1000.\n")
} else {
  cat("\nRows whose sums differ from 1000:\n")
  print(nucleotide_df[different_rows, , drop = FALSE])
  
  cat("\nTheir row sums:\n")
  print(row_totals[different_rows])
}

# Answers to Part 4:
#
# 1. The expected sum is 1000 for each complete kilobase
#    containing only A, C, G, and T.
#
# 2. Yes. Only kb_183 has a sum different from 1000.
#    A = 155, C = 92, G = 80, T = 112.
#    These counts add up to 439.
#
# 3. The sequence contains 182439 bases.
#    There are 182 complete blocks of 1000 bases,
#    followed by a final block containing only 439 bases.
#    This explains why the last row totals 439.
#    Unknown bases such as N could also reduce A/C/G/T
#    totals in general, but that is not the cause here.


# Save the list and data frame for inspection.
saveRDS(kb_counts, "kilobase_counts.rds")

write.csv(
  nucleotide_df,
  "nucleotide_counts.csv",
  row.names = TRUE
)

write.csv(
  data.frame(
    Kilobase = names(kb_counts),
    Sum = row_totals
  ),
  "row_sums.csv",
  row.names = FALSE
)