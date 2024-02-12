function computeScalingFactorB(rateMatrix, baseFreqs) {
    B = 0;
    for (n1 = 0; n1 < Rows(rateMatrix); n1 = n1+1) {
        for (n2 = 0; n2 < Columns(rateMatrix); n2 = n2+1) {
            if (n2 != n1) {
                B = B + baseFreqs[n1]*baseFreqs[n2]*rateMatrix[n1][n2];
            }
        }
    }
    return B;
}

// load the protein alignment including a tree
SetDialogPrompt("Enter path to file containing AA sequences:");
DataSet ds = ReadDataFile(PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter(ds, 1);

// pick one of the AA substitution models
//SelectTemplateModel(filteredData);
#include(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "Jones.mdl");

// input the tree to constrain
IS_TREE_PRESENT_IN_DATA = 0;
ACCEPT_BRANCH_LENGTHS = 1;
ACCEPT_ROOTED_TREES = 1;
SetDialogPrompt("Enter file containing tree:");
fscanf(PROMPT_FOR_FILE, "Raw", treestring);
Tree givenTree = treestring;

// constrain branch lengths 
global scalingB 		= computeScalingFactorB (jonesModel, equalFreqs);
branchNames 	= BranchName (givenTree, -1);
branchLengths	= BranchLength (givenTree, -1);
for (k = 0; k < Columns(branchNames)-1; k=k+1) {
    ExecuteCommands("givenTree." + branchNames[k] + ".t:=" + branchLengths[k] + "/scalingB;");
}

// fit the AA model to the alignment and tree
LikelihoodFunction lf = (filteredData, givenTree);
Optimize (res, lf);

// ancestral reconstruction by joint maximum likelihood
DataSet anc = ReconstructAncestors(lf);
DataSetFilter adf = CreateFilter(anc, 1);

SetDialogPrompt("File name to write ancestral sequences to:");
DATA_FILE_PRINT_FORMAT = 9; // FASTA sequential
fprintf(PROMPT_FOR_FILE, CLEAR_FILE, adf);
