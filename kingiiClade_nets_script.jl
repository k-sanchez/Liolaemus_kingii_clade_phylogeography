#=
# Add packages
Pkg.add("PhyloNetworks") # >= version 0.15
Pkg.add("PhyloPlots") # visualize networks
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("RCall")
Pkg.add("QuartetNetworkGoodnessFit") #  >= version 0.4
Pkg.add("Gadfly")
Pkg.add("JLD") # save Julia session
=#

using Pkg, PhyloNetworks, QuartetNetworkGoodnessFit, PhyloPlots, CSV, DataFrames, JLD, Distributed, Gadfly, RCall; # start a session
Pkg.status("QuartetNetworkGoodnessFit") # check version

cd("PhyloNetworks") # use double quotes!
pwd()

# Add processors
# no needed if you call julia from terminal:
# julia -p 20 script.jl

# addprocs(20)
# @everywhere using PhyloNetworks
# Threads.nthreads() # check n threads

## Input

# CF table
SNPs2CF_file = "SNPs2CF/btw_sp_100quart_boot.csv"
SNPs2CF = readTableCF(SNPs2CF_file)

# guide tree
sppTree_file = "starting_tree.tre";
sppTree = readTopology(sppTree_file);

########################
## Network estimation ##
########################

net0 = snaq!(sppTree, SNPs2CF, hmax = 0, filename = "btw_sp_net0", runs = 10, seed = 233342);
net0 = readSnaqNetwork("btw_sp_net0.out") # optimal net0 can be rooted on silvanae

net1 = snaq!(net0, SNPs2CF, hmax = 1, filename = "btw_sp_net1", runs = 10, seed = 456456)
nets1 = readMultiTopology("btw_sp_net1.out");
net1 = readTopology("(silvanae,((kingii1,#H9:0.0::0.143):0.647,((deseado,(((archef_scol_zully,tristis):0.083,kingii):0.0)#H9:0.0::0.857):0.224,(baguali,escarch_tari):0.321):0.346));")
# optimal net1 cannot be rooted on silvanae, second net has hybrid edge between silvanae and focal group -> I selected nets1[4] (3rd net)

net2 = snaq!(net1, SNPs2CF, hmax = 2, filename = "btw_sp_net2", runs = 10, seed = 898)
nets2 = readMultiTopology("btw_sp_net2.out");
net2 = readTopology("(kingii1,#H9:0.0::0.164,(((baguali,escarch_tari):0.316,(deseado,(((archef_scol_zully,(tristis)#H10:::0.779):0.119,(kingii,#H10:::0.221):0.119):0.0)#H9:0.0::0.836):0.231):0.352,silvanae):0.519);")
# first three nets have 1 hybrid edge -> I selected nets2[5] (4th net)

net3 = snaq!(net2, SNPs2CF, hmax = 3, filename = "btw_sp_net3", runs = 50, seed = 2399898)
nets3 = readMultiTopology("btw_sp_net3.out");
net3 = readTopology("(baguali,(silvanae,(((deseado,(((kingii,#H17:::0.215):0.12,(archef_scol_zully,(tristis)#H17:::0.785):0.12):0.0)#H9:0.0::0.782):0.358,(kingii1,#H9:0.129::0.218):0.453):0.0)#H15:0.0::0.571):0.983,(escarch_tari,#H15:0.233::0.429):0.029);")
# first 48 nets have 1-2 hybrid edges and/or could not be rooted in L. silvanae, selected nets3[49] (4th net)

net4 = snaq!(net3, SNPs2CF, hmax = 4, filename = "btw_sp_net4", runs = 50, seed = 2534542)
nets4 = readMultiTopology("btw_sp_net4.out");
# None of the networks with 4 hybrid edges


# Plot networks

# first and second element of array are the same optimal net (start ploting form [2])
# use this to plot and check the different networks
plot(net1, :R, useEdgeLength = false, style = :majortree, showGamma = true, showEdgeNumber = true);
rootonedge!(net1, 10)

# Final networks (hmax = 0-3)

net0 = readTopology("(baguali,escarch_tari,((kingii1,silvanae):0.2604265565014862,((tristis,archef_scol_zully):0.061082642750672905,(deseado,kingii):0.06783362229427582):0.11458992124437374):0.3770793381525688);")
rootonedge!(net0, 4)
loglik_net0 = 42049.0535425142
net1 = readTopology("(kingii1,#H9:::0.2560620937394804,(((escarch_tari,baguali):0.34993266011230706,((kingii,deseado):0.022923527433140373,(tristis,(archef_scol_zully)#H9:::0.7439379062605196):0.15107483527118604):0.16466830295668464):0.3063237850148288,silvanae):0.39177502982750256);")
rootonedge!(net1, 15)
loglik_net1 = 39091.782999292
net2 = readTopology("(kingii1,#H9:0.0::0.1637996622068637,(((baguali,escarch_tari):0.31607571710366594,(deseado,(((archef_scol_zully,(tristis)#H10:::0.7791274597016159):0.11877243060937308,(kingii,#H10:::0.22087254029838402):0.11877243060937308):0.0)#H9:0.0::0.8362003377931363):0.23075337232457477):0.351574979346867,silvanae):0.5193455157246124);")
rootonedge!(net2, 18)
loglik_net2 = 40878.299818464
net3 = readTopology("(baguali,(silvanae,(((deseado,(((archef_scol_zully,(tristis)#H17:::0.7497135716828062):0.13187968020340377,(kingii,#H17:::0.25028642831719383):0.13187968020340377):0.0)#H9:0.0::0.7671916419836066):0.36824126780781913,(kingii1,#H9:0.00028578619210420554::0.23280835801639344):0.4155930738669343):0.0)#H15:0.0::0.5505420231489968):1.196644613736501,(escarch_tari,#H15:0.3179594776076759::0.44945797685100314):0.14137143877714178);")
rootonedge!(net3, 2)
loglik_net3 = 46789.3171745189

#############################
## Broken stick-heuristics ##
#############################

# best fit for lower values
scores = [loglik_net0, loglik_net1, loglik_net2, loglik_net3]
R"plot"(scores, type = "b", ylab = "network score", xlab = "hmax", col = "blue");


###########
## Plots ##
###########

# rotate edges around some nodes, to plot the networks in the same way
plot(net1, :R, useEdgeLength = false, style = :majortree, showNodeNumber = true); # check
for n in [9, -4, -5] rotate!(net0, n); end # net0
for n in [10, -4, -7] rotate!(net1, n); end # net1
for n in [11, -4, -6] rotate!(net2, n); end # net2
for n in [13] rotate!(net3, n); end # net3

fulltaxnamedict = Dict(
  "silvanae"          =>  "silvanae", # outgroup
  "archef_scol_zully" =>  "archeforus + scolaroi-zullyae",
  "kingii"            =>  "kingii sensu stricto + kingii 2",
  "escarch_tari"      =>  "escarchadosi + tari",
  "kingii1"           =>  "kingii 1",
  "deseado"           =>  "Deseado clade",
  "baguali"           =>  "baguali",
  "tristis"           =>  "tristis",
  )

# Tips with pop/spp names
net0_popNames = deepcopy(net0);
for n in net0_popNames.leaf n.name = fulltaxnamedict[n.name]; end
net1_popNames = deepcopy(net1);
for n in net1_popNames.leaf n.name = fulltaxnamedict[n.name]; end
net2_popNames = deepcopy(net2);
for n in net2_popNames.leaf n.name = fulltaxnamedict[n.name]; end
net3_popNames = deepcopy(net3);
for n in net3_popNames.leaf n.name = fulltaxnamedict[n.name]; end

# net0
R"pdf"("net0Nice.pdf", width = 7, height = 5);
R"par"(mar = [0, 0, 0, 0]);
res = plot(net0_popNames, :R, xlim = [1, 10], tipOffset = 0.1, useEdgeLength = false);
R"dev.off()";

# net1
hi = findfirst([!e.isMajor for e in net1_popNames.edge]) # hi for hybrid index (edge number): 2
R"pdf"("net1Nice.pdf", width = 7, height = 5);
R"par"(mar = [0, 0, 0, 0]);
res = plot(net1_popNames, :R, xlim = [1, 17], style = :majortree, tipOffset = 0.1, useEdgeLength = false, arrowlen = 0.2);
R"text"(res[14][hi, :x] + 0.6, res[14][hi, :y], res[14][hi, :gam], col = "deepskyblue", cex = 0.75);
R"dev.off()";

# net2
hi = findall([!e.isMajor for e in net2_popNames.edge]) # hi for hybrid index
h1 = hi[1] # first minor hybrid edge (edge number): 2
h2 = hi[2] # second minor hybrid edge (edge number): 12
R"pdf"("net2Nice.pdf", width = 7, height = 5);
R"par"(mar = [0, 0, 0, 0]);
res = plot(net2_popNames, :R, xlim = [1.0, 16], style = :majortree, tipOffset = 0.1, useEdgeLength = false);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col = "deepskyblue", length = 0.08, angle = 20);
# R"text"(x = 1, y = 17.5, labels = "no branch lengths", adj = [0,0], cex = 0.8);
x0 = 14.7; x1 = 15.2
R"text"(res[14][h1, :x] + 0.6, res[14][h1, :y], res[14][h1, :gam], col = "deepskyblue", cex = 0.75);
R"text"(res[14][h2, :x] + 0.6, res[14][h2, :y], res[14][h2, :gam], col = "deepskyblue", cex = 0.75);
R"dev.off()";

# net3
# Nice plot, showing topology only (not branch lengths)
hi = findall([!e.isMajor for e in net3_popNames.edge]) # hi for hybrid index
h1 = hi[1] # first minor hybrid edge (edge number): 9
h2 = hi[2] # second minor hybrid edge (edge number): 15
h3 = hi[3] # third minor hybrid edge (edge number): 21
R"pdf"("net3Nice.pdf", width = 7, height = 5);
R"par"(mar = [0, 0, 0, 0]);
res = plot(net3_popNames, :R, xlim = [1.0, 22], tipOffset = 0.1, useEdgeLength = false, style = :majortree);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col = "deepskyblue", length = 0.08, angle = 20);
R"text"(x = 1, y = 17.5, labels = "no branch lengths", adj = [0, 0], cex = 0.8);
x0 = 14.7; x1 = 15.2
R"text"(res[14][h1, :x] + 0.6, res[14][h1, :y], res[14][h1, :gam], col = "deepskyblue", cex = 0.75);
R"text"(res[14][h2, :x] + 0.6, res[14][h2, :y], res[14][h2, :gam], col = "deepskyblue", cex = 0.75);
R"text"(res[14][h3, :x] + 0.6, res[14][h3, :y], res[14][h3, :gam], col = "deepskyblue", cex = 0.75);
R"dev.off()";


###########################
## Goodness of fit tests ##
###########################

## Observed vs expected CFs

using Statistics
cd("PhyloNetworks/goodness_fit/")

topologyMaxQPseudolik!(net3, SNPs2CF); # do this for all networks
df_wide = fittedQuartetCF(SNPs2CF)
df_long = fittedQuartetCF(SNPs2CF, :long)

@rlibrary ggplot2
ggplot(df_long, aes(x = :obsCF, y = :expCF)) + theme_classic() +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "#008080", size = 0.3) + # diagonal line
    geom_point(alpha = 0.5, color = "#008080", position = position_jitter(width = 0.005, height = 0.005)) +
    ylab("quartet CF expected from network") + xlab("quartet CF observed in gene trees") + coord_equal(ratio = 1);
# save with:
ggsave("expCFs_obsVsFittedNet3.svg", scale = 1, width = 6, height = 5);
CSV.write("fitted_CF_net3.csv", df_long)

## Goodness-of-fit test (z statistic)
# nsim = 1000 o mas
# Cai & Ané 2020 -> test statistic: They recommend the use of G (LRT)
# Using simulations I can calculate an empirical p-value as the proportion of simualted z values [6] ⩾ than the observed z [2]

# first, generate a CF table with between-species-only quartets:
# SNPs2CF(), in R, specifying n.quartets = 1
# load CF table as data frame, rename column 'genes' to 'ngenes':
SNPs2CF_DF = DataFrame(CSV.File("../SNPs2CF/btw_sp_1quart.csv"))
rename!(SNPs2CF_DF,:genes => :ngenes)

res0 = quarnetGoFtest!(net0, SNPs2CF_DF, true, seed = 343523, nsim = 100000);
res0[[1, 2, 3]] # p-value, uncorrected z, σ
res0[2]/res0[3] # corrected z
res0[4] # list of outlier p-values (also added to the DF)
res1 = quarnetGoFtest!(net1, SNPs2CF_DF, true, seed = 35354, nsim = 100000);
res1[[1, 2, 3]] # p-value, uncorrected z, σ
res1[2]/res1[3] # corrected z
res2 = quarnetGoFtest!(net2, SNPs2CF_DF, true, seed = 58585, nsim = 100000);
res2[[1, 2, 3]] # p-value, uncorrected z, σ
res2[2]/res2[3] # corrected z
res2[5] # ultrametric network
res3 = quarnetGoFtest!(net3, SNPs2CF_DF, true, seed = 67888698, nsim = 100000);
res3[[1, 2, 3]] # p-value, uncorrected z, σ
res3[2]/res3[3] # corrected z
res3[5] # ultrametric network

# calculate empirical p-value for each network -> proportion of simulated z-values [6] ⩾ observed z-values [2]
using Statistics
# one-sided test: Prob(Z > z)
emp_pvalue0 = mean(sort!(res0[6]) .>= res0[2])
emp_pvalue1 = mean(sort!(res1[6]) .>= res1[2])
emp_pvalue2 = mean(sort!(res2[6]) .>= res2[2])
emp_pvalue3 = mean(sort!(res3[6]) .>= res3[2])

zvalue_observed = res0[2]
zvalue_bootstrap = sort!(res0[6]) # long vector: sorted z-values simulated under the network
pvalue = mean(zvalue_bootstrap .>= zvalue_observed) # one-sided test: Prob(Z > z)

goodFitSummary = DataFrame(network      =   ["net0", "net1", "net2", "net3"],
                           z            =   [res0[2], res1[2], res2[2], res3[2]],
                           sigma        =   [res0[3], res1[3], res2[3], res3[3]],
                           correted_z   =   [res0[2]/res0[3], res1[2]/res1[3], res2[2]/res2[3], res3[2]/res3[3]],
                           pValue       =   [res0[1], res1[1], res2[1], res3[1]],
                           pValueEmp    =   [emp_pvalue0, emp_pvalue1, emp_pvalue2, emp_pvalue3]
                           )
CSV.write("good_fit_QNGF_100Ksims.csv", goodFitSummary)
# CSV.write("SNPs2CF/btw_sp_1quart_pvalue.csv", SNPs2CF_dataFrame) # export DF if desired

## Goodness of fit tests (TICR)

res0 = ticr!(net0, SNPs2CF_DF, true, quartetstat = :maxCF, test = :onesided)
res1 = ticr!(net1, SNPs2CF_DF, true, quartetstat = :maxCF, test = :onesided)
res2 = ticr!(net2, SNPs2CF_DF, true, quartetstat = :maxCF, test = :onesided)
res3 = ticr!(net3, SNPs2CF_DF, true, quartetstat = :maxCF, test = :onesided)

res0[[2, 1]] # z statistic, p-value
res3[3] # count of p-values in each of the four categories

goodFitSummary = DataFrame(network  =   ["net0", "net1", "net2", "net3"],
                           z        =   [res0[2], res1[2], res2[2], res3[2]],
                           pValue   =   [res0[1], res1[1], res2[1], res3[1]]
                           )
CSV.write("good_fit_TICR_multSNPs.csv", goodFitSummary)

outlierPvalues = DataFrame(t1 = SNPs2CF_DF[!, 1],
                           t2 = SNPs2CF_DF[!, 2],
                           t3 = SNPs2CF_DF[!, 3],
                           t4 = SNPs2CF_DF[!, 4],
                           outlierP_net0 = res0[5],
                           outlierP_net1 = res1[5],
                           outlierP_net2 = res2[5],
                           outlierP_net3 = res3[5]
                           )
CSV.write("out_p_values_TICR_multSNPs.csv", outlierPvalues)
