#  MATLAB Code: Gene Expression Data Simulation and Estimation

#  This MATLAB script performs the simulation of genetic data, specifically working with SNP (Single Nucleotide Polymorphism) data and gene expression. The code applies an EM (Expectation-Maximization) algorithm to estimate effect sizes (`alpha`) and variances (`sigma1`, `sigma2`, `sigmau`) using both genotypic and phenotypic data.

#  Key Variables:
# **n1, n2**: Number of samples in the first group (e.g., 205) and second group (e.g., 16,470).
# **m1**: Number of genes (e.g., 440).
# **p1**: Number of SNPs (e.g., 6,769).
# **alpha_initial**: Initial effect size estimates for each gene.
# **k**: SNP loadings, indicating the relative influence of different SNPs on gene expression.
# **wg1_mat, wg2_mat**: Genotypic matrices for two groups, `wg1_mat` for group 1 and `wg2_mat` for group 2.
# **Y_mat**: Gene expression matrix for group 1.
# **newz_mat**: Standardized phenotype data for group 2.

# The script begins by initializing key variables, such as the sample sizes `n1`, `n2`, the number of genes `m1`, the number of SNPs `p1`, and the initial effect sizes `alpha_initial`. SNP loadings `k` are predefined for different genes, representing the number of SNPs that contribute to each gene's expression.

# n1 = 386; n2 = 16470  ; m1 = 507; p1  = 3395; 
# alpha_initial = normrnd(0, 0.01, [m1, 1]);
# old = [0.01,0.01,0.01, alpha_initial'];

# The following is the pathway grouping information:
# k = [repmat(9, 1, 9), repmat(6, 1, 6), repmat(7, 1, 7), repmat(13, 1, 13), ...
#      repmat(9, 1, 9), repmat(18, 1, 18), repmat(16, 1, 16), repmat(7, 1, 7), ...
#      repmat(5, 1, 5), repmat(12, 1, 12), repmat(11, 1, 11), repmat(8, 1, 8), ...
#      repmat(5, 1, 5), repmat(7, 1, 7), repmat(8, 1, 8), repmat(16, 1, 48), ...
#      repmat(9, 1, 9), repmat(7, 1, 7), repmat(5, 1, 5), repmat(6, 1, 6), ...
#      repmat(9, 1, 9), repmat(24, 1, 24), repmat(5, 1, 5), repmat(22, 1, 22), ...
#      repmat(6, 1, 6), repmat(7, 1, 7), repmat(8, 1, 8), repmat(11, 1, 11), ...
#      repmat(5, 1, 10), repmat(11, 1, 11), repmat(5, 1, 5), repmat(13, 1, 13), ...
#      repmat(10, 1, 10), repmat(9, 1, 9), repmat(16, 1, 16), repmat(6, 1, 6), ...
#      repmat(19, 1, 19), repmat(6, 1, 6), repmat(8, 1, 8), repmat(31, 1, 31), ...
#      repmat(9, 1, 18), repmat(5, 1, 5), repmat(6, 1, 6), repmat(5, 1, 5), ...
#      repmat(16, 1, 16), repmat(5, 1, 5)];

#  get initial values
# Y_mat = readtable("y_gene_heart_chu1.csv");
# Y_mat = table2array(Y_mat);

# wg1_mat = readtable('w1_heart_chu1.csv', 'ReadVariableNames', false);
# wg1_array = table2array(wg1_mat);

# non_numeric_columns = [];
# for col = 1:size(wg1_array, 2)
#     if iscell(wg1_array(:, col)) && any(~cellfun(@isnumeric, wg1_array(:, col)))
#         non_numeric_columns = [non_numeric_columns, col];
#     end
# end
# for col = non_numeric_columns
#     wg1_array(:, col) = num2cell(str2double(wg1_array(:, col)));
# end
# wg1_numeric = cell2mat(wg1_array);

# get_mode = @(x) mode(x(~isnan(x)));
# for col = 1:size(wg1_numeric, 2)
#     if any(isnan(wg1_numeric(:, col)))
#         mode_value = get_mode(wg1_numeric(:, col));
#         wg1_numeric(isnan(wg1_numeric(:, col)), col) = mode_value;
#     end
# end

# u = pinv(wg1_numeric' * wg1_numeric) * wg1_numeric' * Y_mat;

# wg2_mat = readtable('w2_heart_chu1.csv', 'ReadVariableNames', false);
# wg2_mat = table2array(wg2_mat);

# X = wg2_mat * u;
# newz_mat = load('z_brain.mat');
# newz_mat  = newz_mat.w1_brain;
# newz_mat = (newz_mat - mean(newz_mat)) ./ std(newz_mat);

# alpha1 = pinv(X' * X) * X' * newz_mat; 

# e1 = Y_mat - wg1_numeric * u;
# e2 = newz_mat - X * alpha1;

# sigma1 = std(e1);
# sigma2 = std(e2);
# sigmau = std(u);
# old = [mean(sigma1), mean(sigma2), mean(sigmau), alpha1'];
# old(1)=0.01;

# a = 0.5;
# lam = 75;
# results = EM_updated(old, a, wg1_numeric, wg2_mat, newz_mat, Y_mat, lam, k, m1, n1, n2, p1);

# The following is the R code.
# In order to save the computation time, the results is saved in the data file called "alphaest_twas1.RData".

# Load necessary data
# Load gene names for the heart tissue dataset
y_gene_heart_chu1_genename <- read_csv("y_gene_heart_chu1_genename.csv")

# Load the pre-saved estimates for TWAS alpha values
load("alphaest_twas1.RData")

# Define key parameters
m1_ch1 = 507    # Number of genes
n1 = 386        # Number of samples in the first group
n2 = 16470      # Number of samples in the second group
p1 = 3395       # Number of SNPs (genetic markers)

# Extract the estimated alpha values (effect sizes) from the loaded data
alphaest_ch1 <- alphaest_twas1[4:(m1_ch1 + 3)]

# Extract the estimated E2 values (model fit residuals) from the loaded data
E2_ch1 <- alphaest_twas1[(m1_ch1 + 5):(m1_ch1 + m1_ch1 + 4)]

# Calculate the mean of the sigma2 estimates
sigma2est = mean(alphaest_twas1[2])

# Define the SNP loadings for each pathway. This vector indicates how many SNPs contribute to each pathway.
k = c(9, 6, 7, 13, 9, 18, 16, 7, 5, 12, 11, 8, 5, 7, 8, 16, 16, 16, 9, 7, 5, 6, 9, 24, 5, 22, 6, 7, 8, 11, 
      5, 5, 11, 5, 13, 10, 9, 16, 6, 19, 6, 8, 31, 9, 9, 5, 6, 5, 16, 5)

# Compute the cumulative sum to identify the SNP block ranges for each pathway
cumulative_sums <- cumsum(k)

# Create starting and ending indices for each pathway block
starts <- c(1, head(cumulative_sums, -1) + 1)
ends <- cumulative_sums

# Create a list of SNP indices for each pathway
group_indices <- mapply(seq, from = starts, to = ends, SIMPLIFY = FALSE)

# Initialize a vector to store p-values for each pathway
p_values_lrt1 = rep(0, length(group_indices))

# Loop over each pathway to calculate p-values using the likelihood ratio test
for (i in 1:length(group_indices)) {
  index = group_indices[[i]]  # SNP indices for the current pathway
  
  # Initialize null hypothesis E2 values (with all pathways as null)
  E2_null = E2_ch1
  
  # Set the E2 values of the current pathway to the null value (mean sigma2 estimate)
  E2_null[index] = alphaest_twas1[2 * m1_ch1 + 5]
  
  # Calculate the log-likelihood for the null model (no effect in the pathway)
  log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum(E2_null) / (length(index) * sigma2est)
  
  # Calculate the log-likelihood for the alternative model (effect in the pathway)
  log_alter = -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum(E2_ch1) / (length(index) * sigma2est)
  
  # Calculate the likelihood ratio statistic and p-value for the current pathway
  lr_statistic = 2 * (log_alter - log_likelihood_z)
  p_values_lrt1[i] = pchisq(lr_statistic, df = length(index) - 1, lower.tail = FALSE)
}

# Initialize a matrix to store the results (pathway names and p-values)
result1 = matrix(NA, 50, 2)

# Define the names of the pathways being tested
pathway_names <- c(
  "Reactome Signaling by PDGF",                                                  
  "Reactome Potassium Channels",                                                 
  "Reactome Signaling by FGFR in disease",                                        
  "KEGG Insulin signaling pathway",                                               
  "KEGG Viral myocarditis",                                                       
  "Reactome Metabolism of carbohydrates",                                         
  "KEGG Alzheimer's disease",                                                     
  "Reactome TRIF mediated TLR3 signaling",                                        
  "Reactome Regulation of Insulin Secretion by Glucagon-like Peptide-",           
  "KEGG Purine metabolism",                                                       
  "Reactome Late Phase of HIV Life Cycle",                                        
  "KEGG Pyruvate metabolism",                                                     
  "KEGG Wnt signaling pathway",                                                   
  "KEGG Oocyte meiosis",                                                          
  "Reactome Unfolded Protein Response",                                           
  "Reactome Phospholipid metabolism",                                             
  "Reactome Innate Immune System",                                                
  "Reactome SLC-mediated transmembrane transport",                                
  "Reactome Mitotic G2-G2/M phases",                                              
  "Reactome MyD88:Mal cascade initiated on plasma membrane",                      
  "Reactome Downstream signal transduction",                                      
  "Reactome G alpha (s) signalling events",                                       
  "KEGG Spliceosome",                                                             
  "Reactome Influenza Viral RNA Transcription and Replication",                   
  "KEGG Vibrio cholerae infection",                                               
  "Reactome Peptide chain elongation",                                            
  "Reactome Recruitment of mitotic centrosome proteins and complexes",            
  "KEGG Propanoate metabolism",                                                   
  "KEGG Peroxisome",                                                              
  "Reactome Cyclin E associated events during G1/S transition",                   
  "KEGG Fc gamma R-mediated phagocytosis",                                        
  "Reactome Lipoprotein metabolism",                                              
  "Reactome Destabilization of mRNA by AUF1 (hnRNP D0)",                          
  "KEGG Type I diabetes mellitus",                                                
  "Reactome Regulation of mitotic cell cycle",                                    
  "Reactome Toll Receptor Cascades",                                              
  "KEGG Glycerophospholipid metabolism",                                          
  "Reactome Regulation of mRNA Stability by Proteins that Bind AU-rich Elements", 
  "Reactome Class B/2 (Secretin family receptors)",                               
  "Reactome Platelet activation, signaling and aggregation",                      
  "Reactome Endosomal Sorting Complex Required For Transport (ESCRT)",            
  "KEGG RNA degradation",                                                         
  "Reactome Cell Cycle",                                                          
  "Reactome Biological oxidations",                                               
  "Reactome Regulation of ornithine decarboxylase (ODC)",                         
  "KEGG Sphingolipid metabolism",                                                 
  "Reactome Transcriptional Regulation of White Adipocyte Differentiation",       
  "KEGG Vasopressin-regulated water reabsorption",                                
  "Reactome Fatty acid, triacylglycerol, and ketone body metabolism",             
  "KEGG Hematopoietic cell lineage"
)

# Store pathway names and corresponding p-values in the result matrix
result1[,1] <- pathway_names 
result1[,2] <- p_values_lrt1



