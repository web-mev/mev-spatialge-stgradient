suppressMessages(suppressWarnings(library('spatialGE')))
suppressMessages(suppressWarnings(library('optparse')))
source('/usr/local/bin/prep_stlist.R')

# args from command line:
args <- commandArgs(TRUE)

# note, -k --kclusters can be integer or character ('dtc') for valid options
option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-c', '--coordinates_file'),
        help='Path to the barcode spatial coordinates input.'
    ),
    make_option(
        c('-s', '--sample_name'),
        help='Sample name'
    ),
    make_option(
        c('-b', '--barcodes'),
        help='barcodes in selected reference cluster'
    ),
    make_option(
        c('-n', '--normalization'),
        help='Normalization method of `log` or `sct`'
    ),
    make_option(
        c('-t', '--topgenes'),
        help='Number of top genes to consider'
    ),
    make_option(
        c('-d', '--distancesummary'),
        help='Distance summary metric to use: min or avg'
    ),
    make_option(
        c('-x', '--xpos_col'),
        help='The column header for the x-position coordinate metadata'
    ),
    make_option(
        c('-y', '--ypos_col'),
        help='The column header for the y-position coordinate metadata'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$coordinates_file)){
    message('Need to provide a count matrix with the -c/--coordinates_file arg.')
    quit(status=1)
}

# transform the name of the normalization scheme:
if (is.null(opt$normalization)){
    message('Need to provide a normalization scheme with the -n/--normalization arg.')
    quit(status=1)
} else if(tolower(opt$normalization) == 'sctransform'){
    norm_scheme <- 'sct'
} else if(tolower(opt$normalization) == 'log'){
    norm_scheme <- 'log'
} else {
    message('We only accept `log` or `SCTransform` for the normalization scheme.')
    quit(status=1)
}

# check distance summary is valid
if (is.null(opt$distancesummary)){
    message('Need to provide a test statistic with the -d/--distancesummary arg.')
    quit(status=1)
} else if(opt$distancesummary == "Minimum"){
    distancesummary <- 'min'
} else if(opt$distancesummary == "Average"){
    distancesummary <- 'avg'
} else {
    message("We only accept `Minimum` or `Average` for the distance summary metric.")
    quit(status=1)
}

tg_as_number <- as.numeric(opt$topgenes)
if (is.na(tg_as_number)){
    message('The -t/--topgenes option must be an integer.')
    quit(status=1)
} else {
    n_topgenes <- as.integer(tg_as_number)
}


# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# call the utility function which will return a list with the necessary items:
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         opt$xpos_col,
                         opt$ypos_col)

# unpack:
spat <- spat_list$spat

# Transform the data
spat <- transform_data(spat, method=norm_scheme)

# Create a pass through of selected barcodes (single cluster) to background pool
ref_barcodes <- make.names(strsplit(opt$barcodes, ',')[[1]])
clusts <- rep(
    2, 
    dim(spat@spatial_meta[[opt$sample_name]])[1]
)
clusts[
    spat@spatial_meta[[opt$sample_name]]$libname %in% ref_barcodes
] <- 1

# Assign passed clusters into object
spat@spatial_meta[[opt$sample_name]]$stclust_pass <- clusts

# Run STgradient
# We can have cores as a param for parallelization "cores = n"
# We can have as input "distsumm = min or avg"
grad_tib <- STgradient(
    spat,
    topgenes = n_topgenes,
    annot = "stclust_pass",
    ref = 1,
    samples = c(opt$sample_name),
    distsumm = distancesummary,
    robust = F,
    cores=2
)

# Write to file
output_filename <- paste(working_dir, 'stgradient_output.tsv', sep='/')
output_table <- grad_tib[[opt$sample_name]]
cn <- colnames(output_table)
# the column names are different depending on the arguments. First remove
# the unnecessary 'sample_name' column and then remove the 'comment' column 
# which is all NA (given the code above)
cn <- cn[cn != 'sample_name']
cn <- cn[!grepl('*_comment$', cn)]
output_table <- output_table[cn]

if(is.null(output_table)){
    message('The output table was empty given the inputs provided. This can happen if the cluster selections and input data are such that a result cannot be calculated.')
    quit(status=1)
}

write.table(
    output_table,
    output_filename,
    sep="\t",
    quote=F,
    row.names=F
)
json_str = paste0('{"STgradient_results":"', output_filename, '"}')
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
