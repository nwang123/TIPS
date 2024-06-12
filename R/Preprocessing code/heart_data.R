heart_geno_gtex <- readRDS("~/Desktop/GTEx_heart_geno_new.Rds")# 838 244276
colnames(heart_geno_gtex) <- gsub("b38.*$", "b38", colnames(heart_geno_gtex))
heart_geneexp <- readRDS("~/Desktop/tpm_heart.Rds") # 2089  435
original_colnames <- colnames(heart_geneexp)
new_names <- gsub("\\.", "-", colnames(heart_geneexp))
new_names <- gsub("^(GTEX-[A-Z0-9]+).*$", "\\1", new_names)
colnames(heart_geneexp) <- new_names
## make gene expression and genotype have the same individuals number n1=386
index1 <- which(heart_geno_gtex$IID %in% colnames(heart_geneexp)[4:435])
heart_geno_gtex <- heart_geno_gtex[index1, ] #386 244276
index2 <- which(colnames(heart_geneexp)[4:435] %in% heart_geno_gtex$IID)
heart_geneexp <- heart_geneexp[,c(1:3,index2+3)] #2089  389

heart_genesnp <- readRDS("~/Desktop/heart_gtex_25kb_genmean_lgmedian_genelist_new.Rds")
gtex_snps_vector <- unlist(lapply(heart_genesnp, function(df) df$gtex_SNP))
gtex_chr <- unlist(lapply(heart_genesnp, function(df) df$avsnp147))
for(i in 1:21) {
  # Dynamically create the column name
  col_name <- paste0("genes", i)
  
  # Use lapply to extract the column from each data frame and unlist to create a vector
  genes_vector <- unlist(lapply(heart_genesnp, function(df) df[[col_name]]))
  
  # Dynamically create the variable name
  var_name <- paste0("genes", i, "_vector")
  
  # Assign the vector to the dynamically created variable name in the global environment
  assign(var_name, genes_vector, envir = .GlobalEnv)
}

genes_vectors <- list()
for(i in 1:21) {
  col_name <- paste0("genes", i)
  genes_vectors[[col_name]] <- unlist(lapply(heart_genesnp, function(df) df[[col_name]]))
}

heart_geno_gtex <- as.matrix(heart_geno_gtex)
heart_geno_gtex  <- rbind(heart_geno_gtex , rep(NA,ncol(heart_geno_gtex ))) #387 244276
for (i in 3:ncol(heart_geno_gtex)) {
  # Find the index of the current SNP in the gtex_snps_vector
  geneindex <- match(colnames(heart_geno_gtex)[i], gtex_snps_vector)
  
  if (!is.na(geneindex)) {  # Ensure geneindex is not NA
    # Initialize a flag to keep track if a match is found
    found <- FALSE
    
    # Iterate through the list of gene vectors
    for (genes_vector in genes_vectors) {
      # Check if the current genes_vector at geneindex is in heart_geneexp$Description
      if (genes_vector[geneindex] %in% heart_geneexp$Description) {
        # If a match is found, update heart_geno_gtex and set found to TRUE
        heart_geno_gtex[387, i] <- genes_vector[geneindex]
        found <- TRUE
        break  # Break the loop since a match is found
      }
    }
    
    # Optionally, handle the case where no match is found across all genes_vectors
    if (!found) {
      # For example, set to NA or perform some other action
      # heart_geno_gtex[387, i] <- NA  # Uncomment or modify as needed
    }
  }
}
sum(is.na(heart_geno_gtex[387,]))
which(is.na(heart_geno_gtex[387,]))[1:10]

matches <- unique(heart_geno_gtex[387,3:ncol(heart_geno_gtex)]) %in% unique(heart_geneexp$Description)  ##2083
matching_values <- unique(heart_geno_gtex[387, 3:ncol(heart_geno_gtex)])[matches]
cols_to_keep <- heart_geno_gtex[387, 3:ncol(heart_geno_gtex)] %in% matching_values
heart_geno_gtex <- heart_geno_gtex[, c(1:2, (3:ncol(heart_geno_gtex))[cols_to_keep])] # 387 226963

matches <- heart_geneexp$Description %in% unique(heart_geno_gtex[387,3:ncol(heart_geno_gtex)]) ##2082
heart_geneexp <- heart_geneexp[matches,] #  2082  389
## step1
n1 <-386; m1<- 2082;p1<- 226963
r_squared_values <- numeric(m1)
adjusted_r_squared_values <- numeric(m1)
# Run the regression for each gene using its cis-SNPs
for (i in 1:m1) {
  index <- which(heart_geno_gtex[387,]==heart_geneexp$Description[i])
  w1 <- heart_geno_gtex[1:(nrow(heart_geno_gtex) - 1), index]  # Extract SNP data for this gene, converting to numeric
  y <- as.vector(unlist(heart_geneexp[i,-c(1:3)]))
  # Fit the linear model
  w1_df <- as.data.frame(w1)
  w1_df[] <- lapply(w1_df, function(x) as.numeric(as.character(x)))
  fit <- lm(y ~ ., data = w1_df)
  # Extract the R-squared value and store it
  r_squared_values[i] <- summary(fit)$r.squared
  adjusted_r_squared_values[i] <- summary(fit)$adj.r.squared
}

# use median r^2 to cut one half
index1 <- which(r_squared_values>median(r_squared_values)) 
new_heart_geneexp <- heart_geneexp[index1, ] # 1041  389
index <- list(rep(NA, nrow(new_heart_geneexp)))
for (i in 1:nrow(new_heart_geneexp)){
  index[[i]] <- which(heart_geno_gtex[387,] == new_heart_geneexp$Description[i])
}
new_heart_geno_gtex <- heart_geno_gtex[, c(1:3, unlist(index))] #387 192866


## step 3 use lm model to pick top 10 snp and then use again to check the r^2 and adjusted r^2
step_3 <- function(p){
  top_10_snps_per_gene <- list()
  for (i in 1:nrow(new_heart_geneexp)) {
    index <- which(new_heart_geno_gtex[387,]==new_heart_geneexp$Description[i])
    w1 <- new_heart_geno_gtex[1:(nrow(new_heart_geno_gtex) - 1), index]  # Extract SNP data for this gene, converting to numeric
    y <- as.vector(unlist(new_heart_geneexp[i,-c(1:3)]))
    # Fit the linear model
    
    w1_df <- as.data.frame(w1)
    w1_df[] <- lapply(w1_df, function(x) as.numeric(as.character(x)))
    fit <- lm(y ~ ., data = w1_df)
    # Identify top 10 SNPs based on lowest p-values
    # Extract coefficients summary excluding the intercept
    coef_summary <- summary(fit)$coefficients[-1,]  # Excludes intercept row
    if (is.vector(coef_summary)){
      top_10_snps_per_gene[[i]] <- NA
    }else{
      # Order SNPs by p-value and select top 10
      top_10_indices <- order(coef_summary[, 4])[1:p]
      top_10_indices <- top_10_indices[!is.na(top_10_indices)]
      # Get the SNP names for the top 10 indices
      top_10_snps_per_gene[[i]] <- rownames(coef_summary)[top_10_indices]
    }
    
  }
  
  ## new genotype matrix 
  index <- list(rep(NA, nrow(new_heart_geneexp)))
  for (i in 1:nrow(new_heart_geneexp)){
    index[[i]] <- which(colnames(new_heart_geno_gtex) %in% top_10_snps_per_gene[[i]])
  }
  new_heart_geno_gtex2 <- new_heart_geno_gtex[, c(1:3, unlist(index))]
  
  r_squared_values_new <- numeric(length(top_10_snps_per_gene))
  adjusted_r_squared_values_new <- numeric(length(top_10_snps_per_gene))
  for (i in 1:length(top_10_snps_per_gene)) {
    index <-  which(colnames(new_heart_geno_gtex2) %in% top_10_snps_per_gene[[i]])
    w1 <- new_heart_geno_gtex2[1:(nrow(new_heart_geno_gtex2) - 1), index]  # Extract SNP data for this gene, converting to numeric
    y <- as.vector(unlist(new_heart_geneexp[i,-c(1:3)]))
    # Fit the linear model
    if (length(index)==0){
      w1 <- rep(0,length(y))
    }
    w1_df <- as.data.frame(w1)
    w1_df[] <- lapply(w1_df, function(x) as.numeric(as.character(x)))
    fit <- lm(y ~ ., data = w1_df)
    r_squared_values_new[i] <- summary(fit)$r.squared
    adjusted_r_squared_values_new[i] <- summary(fit)$adj.r.squared
  }
  return(list(r_squared=r_squared_values_new, adjusted_r_squared = adjusted_r_squared_values_new))
}

top10 <- step_3(10)
top20 <- step_3(20)
top30 <- step_3(30)

r_squared_10 <- top10$r_squared;adjusted_r_squared_10 <- top10$adjusted_r_squared
r_squared_20 <- top20$r_squared;adjusted_r_squared_20 <- top20$adjusted_r_squared
r_squared_30 <- top30$r_squared;adjusted_r_squared_30 <- top30$adjusted_r_squared

true <- r_squared_values[which(r_squared_values>median(r_squared_values))]
true_adj <- adjusted_r_squared_values[which(r_squared_values>median(r_squared_values))]

# Assuming 'r_squared_values' and 'adjusted_r_squared_values' are your true values

# Differences in R^2
diff_r_squared_10 <- true - r_squared_10
diff_r_squared_20 <- true - r_squared_20
diff_r_squared_30 <- true - r_squared_30

# Differences in Adjusted R^2
diff_adjusted_r_squared_10 <- true_adj - adjusted_r_squared_10
diff_adjusted_r_squared_20 <- true_adj - adjusted_r_squared_20
diff_adjusted_r_squared_30 <- true_adj - adjusted_r_squared_30

# Proportions of R^2
prop_r_squared_10 <- r_squared_10 / true
prop_r_squared_20 <- r_squared_20 / true
prop_r_squared_30 <- r_squared_30 / true



summarize_metrics <- function(values) {
  summary_values <- summary(values)
  summary_values[c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")]
}

# Summarize differences in R^2
cat("Differences in R^2\n")
cat("Top 10:\n"); print(summarize_metrics(diff_r_squared_10))
cat("Top 20:\n"); print(summarize_metrics(diff_r_squared_20))
cat("Top 30:\n"); print(summarize_metrics(diff_r_squared_30))

cat("\nDifferences in Adjusted R^2\n")
cat("Top 10:\n"); print(summarize_metrics(diff_adjusted_r_squared_10))
cat("Top 20:\n"); print(summarize_metrics(diff_adjusted_r_squared_20))
cat("Top 30:\n"); print(summarize_metrics(diff_adjusted_r_squared_30))

# Summarize proportions of R^2
cat("\nProportions of R^2\n")
cat("Top 10:\n"); print(summarize_metrics(prop_r_squared_10))
cat("Top 20:\n"); print(summarize_metrics(prop_r_squared_20))
cat("Top 30:\n"); print(summarize_metrics(prop_r_squared_30))

length(which(new_heart_geneexp$Description %in% y_gene$Description))



## total 1041 genes
high_r2_genes <- which(prop_r_squared_20 > 0.5)
current_r_squared <- matrix(NA,20,1041)
for(i in 1:20) {  # Assuming you start with top 20 SNPs and go down to 1
  current_r_squared[i,] <- step_3(i)$r_squared
}

min_snps_required <- rep(NA, 1041)
for(gene in high_r2_genes) {
  for(num_snps in 1:20) {  # Assuming you start with top 20 SNPs and go down to 1
    current_r_squared1 <- current_r_squared[num_snps,gene] 
    if(current_r_squared1/true[gene] >= 0.5) {
      min_snps_required[gene] <- num_snps
      break  # Exit the loop once the condition is met
    }
  }
}

min_snps_required[is.na(min_snps_required)] <- 20
sum(min_snps_required)

topsnp <- function(max_top){
  # Initialize a list to hold the top SNPs information for each gene and each top size
  top_snps_per_gene_and_size <- list()
  
  for (i in 1:nrow(new_heart_geneexp)) {
    # Find the SNPs corresponding to the current gene
    index <- which(new_heart_geno_gtex[387,] == new_heart_geneexp$Description[i])
    if (length(index) == 0) {
      print(paste("No SNPs found for gene:", new_heart_geneexp$Description[i]))
      next # Skip if no SNPs found for this gene
    }
    
    w1 <- new_heart_geno_gtex[1:(nrow(new_heart_geno_gtex) - 1), index]  # Extract SNP data for this gene
    y <- as.vector(unlist(new_heart_geneexp[i, -c(1:3)]))
    
    w1_df <- as.data.frame(w1)
    colnames(w1_df) <- make.unique(names(w1_df))
    w1_df[] <- lapply(w1_df, as.numeric)
    fit <- lm(y ~ ., data = w1_df)
    
    # Initialize a list to store top SNPs for different top sizes for the current gene
    top_snps_per_size <- list()
    
    # Extract coefficients summary excluding the intercept
    coef_summary <- summary(fit)$coefficients[-1,]  # Excludes intercept row
    if (!is.matrix(coef_summary)) {
      print(paste("No significant SNPs for gene:", new_heart_geneexp$Description[i]))
      next # Skip if coef_summary is not a matrix (e.g., no SNPs were significant)
    } else {
      # Loop through 1 to max_top to get top N SNPs
      for(p in 1:max_top) {
        if (ncol(coef_summary) < 4 || nrow(coef_summary) < p) {
          top_snps_per_size[[as.character(p)]] <- NA # Handle cases with fewer SNPs than p
        } else {
          top_10_indices <- order(coef_summary[, 4], decreasing = FALSE)[1:p]
          top_10_indices <- top_10_indices[!is.na(top_10_indices)]
          top_snps_per_size[[as.character(p)]] <- rownames(coef_summary)[top_10_indices]
        }
      }
    }
    # Store the top SNPs for this gene in the main list
    top_snps_per_gene_and_size[[new_heart_geneexp$Description[i]]] <- top_snps_per_size
  }
  return(top_snps_per_gene_and_size)
}

# Now call the function to get top 1 to 20 SNPs for each gene
top_snps_list <- topsnp(20)
new_heart_geneexp <- new_heart_geneexp[!grepl("EXOSC6", new_heart_geneexp$Description), ]
index <- list(rep(NA, nrow(new_heart_geneexp)))
for (i in 1:nrow(new_heart_geneexp)){
  index[[i]] <- which(colnames(new_heart_geno_gtex) %in% top_snps_list[[i]][[min_snps_required[i]]])
}
new_heart_geno_gtex_final <- new_heart_geno_gtex[, c(1:3, unlist(index))] # 387 17574
index <- vector("list", nrow(new_heart_geneexp))
# Loop through each row of new_heart_geneexp
for (i in 1:nrow(new_heart_geneexp)) {
  # Find the index where the condition is true
  idx <- which(new_heart_geno_gtex_final[387,-c(1:2)] == new_heart_geneexp$Description[i])
  
  # Check if idx is NA
  if (length(idx) == 0) {
    index[[i]] <- NA  # If idx is NA, assign NA to the list element
  } else {
    index[[i]] <- idx  # Otherwise, assign idx to the list element
  }
}
na_indices <- vector()
# Loop through the list and check for NA values
for (i in 1:length(index)) {
  if (any(is.na(index[[i]]))) {
    na_indices <- c(na_indices, i)  # Append the index to the na_indices vector
  }
}
new_heart_geneexp <- new_heart_geneexp[-na_indices,]
new_heart_geno_gtex_final <- new_heart_geno_gtex_final[,c(1,2,unlist(index))]#387 17573
# Apply the function to replace NA values with the mode
for (i in 3:ncol(new_heart_geno_gtex_final)) {
  if (any(is.na(new_heart_geno_gtex_final[, i]))) {
    mode_value <- get_mode(new_heart_geno_gtex_final[, i])
    new_heart_geno_gtex_final[is.na(new_heart_geno_gtex_final[, i]), i] <- mode_value
  }
}

# Check the dimensions
dim(new_heart_geno_gtex_final)
length_heart <-  sapply(index, length)

write.csv(length_heart, "length_heart.csv", row.names = FALSE, col.names = FALSE)
write.csv(t(new_heart_geneexp[,-c(1:3)]), "y_heart_unique.csv", row.names=FALSE,col.names = FALSE)
write.table(new_heart_geno_gtex_final[-387, -c(1:2)], "w1_heart_unique.csv", row.names=FALSE,col.names = FALSE)

### pathway
load("~/Downloads/Reactome.RData")
load("~/Downloads/KEGG_MSigDB.RData")
load("~/Downloads/biocarta.RData")
reactomeDB <- c(reactomeDB,KEGG_MSigDB,biocartaDB)
matched_genes <- lapply(reactomeDB, function(pathway_genes) {
  indices <- which(pathway_genes %in% new_heart_geneexp$Description )
  list(indices = indices, genes = pathway_genes[indices])
})
names(matched_genes) <- names(reactomeDB)
# combined_genes <- unlist(lapply(matched_genes, function(x) x$genes))
gene_lengths <- sapply(matched_genes, function(pathway) length(pathway$genes))
filtered_gene_lengths <- gene_lengths[gene_lengths >= 5]
filtered_matched_genes <- matched_genes[names(filtered_gene_lengths)]
all_genes <- unlist(sapply(filtered_matched_genes, function(pathway) pathway$genes)) # 2968
names(all_genes) <- gsub("([0-9]+)$", "", names(all_genes))
# Find genes that appear more than once
unique_genes <- unique(all_genes)
element_lengths <- as.vector(sapply(filtered_matched_genes, function(x) length(x$genes)))
k_vec <- c(unlist(sapply(element_lengths, function(x) rep((x), x))),
           rep(1,nrow(new_heart_geneexp)-length(unique_genes)))  # 2866:430 unqiue in pathway, 2256 total in pathway, 610 scattered
## y order with pathway
matched_indices <- rep(0,length(all_genes))
for (i in 1:length(all_genes)){
  matched_indices[i] <- which(new_heart_geneexp$Description %in% all_genes[i])
}
all_indices <- seq_along(new_heart_geneexp$Description)
ordered_indices <- c(matched_indices, setdiff(all_indices, matched_indices))

# Reorder the 'y_filtered$Description' according to the 'ordered_indices'
y_filtered_ordered <- new_heart_geneexp[ordered_indices,] #  2866  389

## new genotype matrix 
index <- list(rep(NA, nrow(y_filtered_ordered)))
for (i in 1:nrow(y_filtered_ordered)){
  index[[i]] <- which(new_heart_geno_gtex_final[387,] == y_filtered_ordered$Description[i])
}
geno_heart <- new_heart_geno_gtex_final[, c(1:2, unlist(index))] # 387 48970
index <- list(rep(NA, nrow(y_filtered_ordered)))
for (i in 1:nrow(y_filtered_ordered)){
  index[[i]] <- which(geno_heart[387,] == y_filtered_ordered$Description[i])
}

gene_geno_index <- (as.vector(geno_heart[387,3:ncol(geno_heart)]))
gene_geno_index <- as.character(gene_geno_index)
rle_result <- rle(gene_geno_index)
index_gene_geno <- rle_result$lengths
write.csv(index_gene_geno, "index_gene_geno.csv", row.names = FALSE, col.names = FALSE)
write.csv(k_vec, "k_vec.csv", row.names = FALSE, col.names = FALSE)


write.csv(w1_heart, "w1_heart.csv", row.names=FALSE,col.names = FALSE)
# Define the get_mode function
get_mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

# Create a copy of geno_heart excluding the 387th row to modify
w1_heart <- geno_heart[-387,]

# Apply the function to replace NA values with the mode
for (i in 3:ncol(w1_heart)) {
  if (any(is.na(w1_heart[, i]))) {
    mode_value <- get_mode(w1_heart[, i])
    w1_heart[is.na(w1_heart[, i]), i] <- mode_value
  }
}

# Add back the first two columns
w1_heart <- cbind(geno_heart[-387, 1:2], w1_heart[, 3:ncol(w1_heart)])

# Check the dimensions
dim(w1_heart)



write.table(w1_heart, 
            "w1_brain.csv", 
            sep = ",",  # Set comma as the field separator
            row.names = FALSE, 
            col.names = FALSE) 
w1_heart <- read.csv('w1_brain.csv', header = FALSE)
write.table(t(y_filtered_ordered[,-c(1:3)]),"y_gene_heart.csv", sep = ",",  row.names = FALSE) 

y <- read.csv("y_gene_heart.csv")
indexgeno <- rep(NA, 17571)
# Loop through each column in geno_heart (excluding the first two columns)
for (i in 1:17571) {
  match_index <- which(colnames(new_heart_geno_gtex_final)[-(1:2)][i] == colnames(heart_geno_gtex)[-(1:2)])
  if (length(match_index) == 0) {
    # If no match is found, set indexgeno[i] to NA or some other placeholder
    indexgeno[i] <- NA
  } else {
    # If a match is found, store the index
    indexgeno[i] <- match_index
  }
}
indexgeno <-  indexgeno[-which(is.na(indexgeno))]
write.csv(indexgeno, "indexgeno.csv", row.names = FALSE, col.names = FALSE)


##w2 
heart_geno_ukb <- readRDS("UKB_heart_geno_new.Rds") # 16470 244276
indexgeno <- as.vector(read.table("indexgeno.csv",header=TRUE)$x)
heart_geno_ukb <- heart_geno_ukb[c(1,2,indexgeno+2)] #16470 48970
# replace na with mode
get_mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

# Applying the mode replacement to each column
heart_geno_ukb <- data.frame(lapply(heart_geno_ukb, function(x) {
  if(any(is.na(x))) {  # Check if there is any NA in the column
    mode_value <- get_mode(x)  # Get mode of the column
    x[is.na(x)] <- mode_value  # Replace NA with the mode
  }
  return(x)
}))
write.csv(heart_geno_ukb[-c(1,2)], "w2_heart.csv", row.names=FALSE)  #16470 48970
write.csv(heart_geno_ukb[-c(1,2)], "w2_heart_unique.csv", row.names=FALSE)  #16470 17571

heart_geno_ukb <- read.csv("w2_heart.csv",header=TRUE) # 16470 48968
heart_pheno <- readRDS("~/Desktop/SBP.Rds") # 16470     2
write.csv(heart_pheno[,-1], "z_heart.csv", row.names=FALSE)



### scatted gene 5 chunks
# Define the start and end rows
start_row <- 2257
end_row <- 2866
# 2866:430 unqiue in pathway, 2256 total in pathway, 610 scattered
# Calculate the number of rows
n_rows <- end_row - start_row + 1

# Generate shuffled indices from start_row to end_row
shuffled_indices <- sample(start_row:end_row)

# Split indices into 5 random chunks
chunks <- split(shuffled_indices, cut(shuffled_indices,6, labels = FALSE))

# You can access each chunk using chunks[[1]], chunks[[2]], etc.
chunked_data <- lapply(chunks, function(chunk) y_filtered_ordered[chunk, ])

get_genotype_data <- function(chunk) {
  # Collect unique descriptions from the chunk
  descriptions <- (chunk$Description)
  # Initialize an empty list to store indices
  index <- list(rep(NA, length(descriptions)))
  # Loop through unique descriptions to find matching indices in the genotype matrix
  for (i in 1:length(descriptions)){
    index[[i]] <- which(new_heart_geno_gtex_final[387,] == descriptions[i])
  }
  
  # Extract these columns to form the genotype data for the chunk
  genotype_chunk <- new_heart_geno_gtex_final[, c(1:2,unlist(index))]
  return(genotype_chunk)
}

# Apply this function to each chunk
geno_heart_chunks1 <- get_genotype_data(chunked_data[[1]])
geno_heart_chunks2 <- get_genotype_data(chunked_data[[2]])
geno_heart_chunks3 <- get_genotype_data(chunked_data[[3]])
geno_heart_chunks4 <- get_genotype_data(chunked_data[[4]])
geno_heart_chunks5 <- get_genotype_data(chunked_data[[5]])
geno_heart_chunks6 <- get_genotype_data(chunked_data[[6]])
geno_heart_scatter <- get_genotype_data(y_filtered_ordered[2257: 2866,])

### all pathway
group_indices <- 1:length(element_lengths)

# Classify groups by size
group_levels <- cut(element_lengths, breaks = c(0, 6, 11, Inf), right = FALSE,
                    labels = c("5-6", "7-11", "12+"))

# Split groups by level
levels_list <- split(group_indices, group_levels)

# Shuffle and chunk
set.seed(42)  # For reproducibility
chunked_groups <- lapply(levels_list, function(indices) {
  shuffled <- sample(indices)
  split(shuffled, ceiling(seq_along(shuffled)/length(shuffled)*6))
})

# Combine chunks from each level
final_chunks <- vector("list", 6)
for (i in 1:6) {
  final_chunks[[i]] <- unlist(lapply(chunked_groups, `[[`, i))
}

# Optionally shuffle each final chunk
final_chunks <- lapply(final_chunks, sample)

# Print indices and corresponding values from element_lengths for each chunk
final_chunks <- lapply(final_chunks, function(chunk) {
  list(
    indices = chunk,
    values = element_lengths[chunk]
  )
})
gene_pathway <- rep(NA,247)
for (i in 1:length(element_lengths)){
  gene_pathway[i] <- names(all_genes)[cumsum(element_lengths)[i]-1]
}

# Initialize an empty list for storing grouped genes
gene_groups <- vector("list", length = length(element_lengths))

# Fill each group based on element_lengths
for (i in 1:length(element_lengths)) {
  if (i==1){
    start_index <- 1
  }else{
    start_index <- cumsum(element_lengths)[i-1]+1 
  }
  end_index <- cumsum(element_lengths)[i] 
  gene_groups[[i]] <- all_genes[start_index:end_index]
}
# Fetching gene information from the gene_groups
genes_chunk1 <- unlist(gene_groups[as.numeric(final_chunks[[1]]$indices)])
genes_chunk2 <- unlist(gene_groups[as.numeric(final_chunks[[2]]$indices)])
genes_chunk3 <- unlist(gene_groups[as.numeric(final_chunks[[3]]$indices)])
genes_chunk4 <- unlist(gene_groups[as.numeric(final_chunks[[4]]$indices)])
genes_chunk5 <- unlist(gene_groups[as.numeric(final_chunks[[5]]$indices)])
genes_chunk6 <- unlist(gene_groups[as.numeric(final_chunks[[6]]$indices)])

# Calculate the size of pathway for each chunk
group_lengths1 <- sapply(gene_groups[as.numeric(final_chunks[[1]]$indices)], length)
group_lengths2 <- sapply(gene_groups[as.numeric(final_chunks[[2]]$indices)], length)
group_lengths3 <- sapply(gene_groups[as.numeric(final_chunks[[3]]$indices)], length)
group_lengths4 <- sapply(gene_groups[as.numeric(final_chunks[[4]]$indices)], length)
group_lengths5 <- sapply(gene_groups[as.numeric(final_chunks[[5]]$indices)], length)
group_lengths6 <- sapply(gene_groups[as.numeric(final_chunks[[6]]$indices)], length)

pathway_gene <- function(all_genes){
  matched_indices <- rep(0,length(all_genes))
  for (i in 1:length(all_genes)){
    matched_indices[i] <- which(new_heart_geneexp$Description %in% all_genes[i])
  }
  all_indices <- seq_along(new_heart_geneexp$Description)
  ordered_indices <- c(matched_indices, setdiff(all_indices, matched_indices))
  
  # Reorder the 'y_filtered$Description' according to the 'ordered_indices'
  y_filtered_ordered <- new_heart_geneexp[matched_indices,]
  return(y_filtered_ordered)
}
y_filtered_ordered1 <- pathway_gene(genes_chunk1)
geno_heart_chunk1 <- get_genotype_data(y_filtered_ordered1)
y_filtered_ordered2 <- pathway_gene(genes_chunk2)
geno_heart_chunk2 <- get_genotype_data(y_filtered_ordered2)
y_filtered_ordered3 <- pathway_gene(genes_chunk3)
geno_heart_chunk3 <- get_genotype_data(y_filtered_ordered3)
y_filtered_ordered4 <- pathway_gene(genes_chunk4)
geno_heart_chunk4 <- get_genotype_data(y_filtered_ordered4)
y_filtered_ordered5 <- pathway_gene(genes_chunk5)
geno_heart_chunk5 <- get_genotype_data(y_filtered_ordered5)
y_filtered_ordered6 <- pathway_gene(genes_chunk6)
geno_heart_chunk6 <- get_genotype_data(y_filtered_ordered6)

y_heart_chu1 <- rbind(y_filtered_ordered1[,],chunked_data[[1]][,])
y_heart_chu2 <- rbind(y_filtered_ordered2[,],chunked_data[[2]][,])
y_heart_chu3 <- rbind(y_filtered_ordered3[,],chunked_data[[3]][,])
y_heart_chu4 <- rbind(y_filtered_ordered4[,],chunked_data[[4]][,])
y_heart_chu5 <- rbind(y_filtered_ordered5[,],chunked_data[[5]][,])
y_heart_chu6 <- rbind(y_filtered_ordered6[,],chunked_data[[6]][,])

write.table(t(y_heart_chu1[,-c(1:3)]),"y_gene_heart_chu1.csv", sep = ",",  row.names = FALSE) 
write.table(t(y_heart_chu2[,-c(1:3)]),"y_gene_heart_chu2.csv", sep = ",",  row.names = FALSE) 
write.table(t(y_heart_chu3[,-c(1:3)]),"y_gene_heart_chu3.csv", sep = ",",  row.names = FALSE) 
write.table(t(y_heart_chu4[,-c(1:3)]),"y_gene_heart_chu4.csv", sep = ",",  row.names = FALSE) 
write.table(t(y_heart_chu5[,-c(1:3)]),"y_gene_heart_chu5.csv", sep = ",",  row.names = FALSE) 
write.table(t(y_heart_chu6[,-c(1:3)]),"y_gene_heart_chu6.csv", sep = ",",  row.names = FALSE) 


geno_heart_chunk1 <- cbind(geno_heart_chunk1,geno_heart_chunks1)
geno_heart_chunk2 <- cbind(geno_heart_chunk2,geno_heart_chunks2)
geno_heart_chunk3 <- cbind(geno_heart_chunk3,geno_heart_chunks3)
geno_heart_chunk4 <- cbind(geno_heart_chunk4,geno_heart_chunks4)
geno_heart_chunk5 <- cbind(geno_heart_chunk5,geno_heart_chunks5)
geno_heart_chunk6 <- cbind(geno_heart_chunk6,geno_heart_chunks6)
geno_heart_chunk1 <- geno_heart_chunk1[, !colnames(geno_heart_chunk1) %in% c("FID", "IID")]
geno_heart_chunk2 <- geno_heart_chunk2[, !colnames(geno_heart_chunk2) %in% c("FID", "IID")]
geno_heart_chunk3 <- geno_heart_chunk3[, !colnames(geno_heart_chunk3) %in% c("FID", "IID")]
geno_heart_chunk4 <- geno_heart_chunk4[, !colnames(geno_heart_chunk4) %in% c("FID", "IID")]
geno_heart_chunk5 <- geno_heart_chunk5[, !colnames(geno_heart_chunk5) %in% c("FID", "IID")]
geno_heart_chunk6 <- geno_heart_chunk6[, !colnames(geno_heart_chunk6) %in% c("FID", "IID")]
write.table(geno_heart_chunk1[-387,], "w1_heart_chu1.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno_heart_chunk2[-387,], "w1_heart_chu2.csv", sep = ",",row.names = FALSE, col.names =  FALSE) 
write.table(geno_heart_chunk3[-387,], "w1_heart_chu3.csv", sep = ",",row.names = FALSE, col.names =  FALSE) 
write.table(geno_heart_chunk4[-387,], "w1_heart_chu4.csv", sep = ",",row.names = FALSE, col.names =  FALSE) 
write.table(geno_heart_chunk5[-387,], "w1_heart_chu5.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno_heart_chunk6[-387,], "w1_heart_chu6.csv", sep = ",",row.names = FALSE, col.names =  FALSE) 

geno_heart_chunk11 <-  data.frame(lapply(geno_heart_chunk1, function(x) {
    if(any(is.na(x))) {  # Check if there is any NA in the column
      mode_value <- get_mode(x)  # Get mode of the column
      x[is.na(x)] <- mode_value  # Replace NA with the mode
    }
    return(x)
  }))


index_w1_w2 <- function(geno_heart_chunk1, heart_geno_gtex) {
  # Create a vector to store the indices
  indexgeno <- rep(NA, ncol(geno_heart_chunk1) )  # Assuming the first 2 cols are not part of the matching
  
  # Loop through the columns of geno_brain_chunk1, excluding the first two columns
  for (i in 1:ncol(geno_heart_chunk1)) {
    # Find the matching index in the brain_geno_gtex (excluding the first two columns)
    match_index <- which(colnames(geno_heart_chunk1)[i] == colnames(heart_geno_gtex))
    
    # Check if match_index is not empty
    if (length(match_index) > 0) {
      indexgeno[i ] <- match_index  # Adjust index to fill from the start of indexgeno
    }  # If match_index is empty, indexgeno[i - 2] remains NA
  }
  
  return(indexgeno)
}

index_w1w21 <- index_w1_w2(geno_heart_chunk1, heart_geno_gtex)
valid_indices1 <- which(cols_to_keep==TRUE)[index_w1w21[!is.na(index_w1w21)]]
write.csv(valid_indices1, "w1w2index1.csv", row.names = FALSE, col.names = FALSE)
index_w1w22 <- index_w1_w2(geno_heart_chunk2, heart_geno_gtex)
valid_indices2 <- which(cols_to_keep==TRUE)[index_w1w22[!is.na(index_w1w22)]]
write.csv(valid_indices2, "w1w2index2.csv", row.names = FALSE, col.names = FALSE)
index_w1w23 <- index_w1_w2(geno_heart_chunk3, heart_geno_gtex)
valid_indices3 <- which(cols_to_keep==TRUE)[index_w1w23[!is.na(index_w1w23)]]
write.csv(valid_indices3, "w1w2index3.csv", row.names = FALSE, col.names = FALSE)
index_w1w24 <- index_w1_w2(geno_heart_chunk4, heart_geno_gtex)
valid_indices4 <- which(cols_to_keep==TRUE)[index_w1w24[!is.na(index_w1w24)]]
write.csv(valid_indices4, "w1w2index4.csv", row.names = FALSE, col.names = FALSE)
index_w1w25 <- index_w1_w2(geno_heart_chunk5, heart_geno_gtex)
valid_indices5 <- which(cols_to_keep==TRUE)[index_w1w25[!is.na(index_w1w25)]]
write.csv(valid_indices5, "w1w2index5.csv", row.names = FALSE, col.names = FALSE)
index_w1w26 <- index_w1_w2(geno_heart_chunk6, heart_geno_gtex)
valid_indices6 <- which(cols_to_keep==TRUE)[index_w1w26[!is.na(index_w1w26)]]
write.csv(valid_indices6, "w1w2index6.csv", row.names = FALSE, col.names = FALSE)

valid_indices1 <- as.vector(read.table("w1w2index1.csv",header=TRUE)$x) 
geno2_heart_chunk1 <- heart_geno_ukb[,valid_indices1]
valid_indices2 <- as.vector(read.table("w1w2index2.csv",header=TRUE)$x) 
geno2_heart_chunk2 <- heart_geno_ukb[,valid_indices2]
valid_indices3 <- as.vector(read.table("w1w2index3.csv",header=TRUE)$x) 
geno2_heart_chunk3 <- heart_geno_ukb[,valid_indices3]
valid_indices4 <- as.vector(read.table("w1w2index4.csv",header=TRUE)$x) 
geno2_heart_chunk4 <- heart_geno_ukb[,valid_indices4]
valid_indices5 <- as.vector(read.table("w1w2index5.csv",header=TRUE)$x) 
geno2_heart_chunk5 <- heart_geno_ukb[,valid_indices5]
valid_indices6 <- as.vector(read.table("w1w2index6.csv",header=TRUE)$x) 
geno2_heart_chunk6 <- heart_geno_ukb[,valid_indices6]

#16470 48968
write.table(geno2_heart_chunk1, "w2_heart_chu1.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno2_heart_chunk2, "w2_heart_chu2.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno2_heart_chunk3, "w2_heart_chu3.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno2_heart_chunk4, "w2_heart_chu4.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno2_heart_chunk5, "w2_heart_chu5.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
write.table(geno2_heart_chunk6, "w2_heart_chu6.csv", sep = ",",row.names = FALSE, col.names = FALSE) 
