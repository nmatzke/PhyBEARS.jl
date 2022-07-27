

wd = "/GitHub/BioGeoJulia.jl/notes/work_precision_fig/"
setwd(wd)

df_fns = c("workprecision_020areas_040states.Rdata",
"workprecision_030areas_080states.Rdata",
"workprecision_040areas_0160states.Rdata",
"workprecision_050areas_0320states.Rdata",
"workprecision_060areas_0640states.Rdata",
"workprecision_070areas_01280states.Rdata",
"workprecision_080areas_02560states.Rdata",
"workprecision_090areas_05120states.Rdata",
"workprecision_0100areas_010240states.Rdata",
"workprecision_0110areas_020480states.Rdata")



k = 1

for (k in 1:length(df_fns))
	{
	# Loads to df
	load(df_fns[k])

	if (k == 1)
		{
		solver_names = names(df)
		dims = c(dim(df), length(df_fns))
		big_cube = array(data=NA, dim=dims)
		}

	for (j in 1:length(solver_names))
		{
		if (solver_names[j] %in% names(df))
			{
			tmp_col = df[,solver_names[j]]
			} else {
			tmp_col = rep(NaN, times=dims[1])
			}
		
		big_cube[,j,k] = tmp_col
		} # end j

	} # end k


#######################################################
# Plot work-precision diagram
#######################################################
pdffn = "work_precision_figs_v1.pdf"
pdf(file=pdffn, width=6, height=6)


for (k in 1:dim(big_cube)[3])
	{
	tmpmat = big_cube[,,k]
	abstol = tmpmat[,1]
	colors = rainbow(n=ncol(tmpmat)-1)
	ticklength = 0.5
	
	tmpfun <- function(tmpcol)
		{
		all(is.nan(tmpcol))
		}
	keep_cols_TF = apply(X=tmpmat[,2:ncol(tmpmat)], MARGIN=2, FUN=tmpfun) == FALSE
	
	plot(x=abstol, y=seq(0.001, 100000, length.out=length(abstol)), xlim=c(max(abstol), min(abstol)), pch=".", col="white", xlab="", ylab="", log="xy", yaxt="n")

	axis(2, at=c(0.001, 0.01, 0.1, 1, 10, 100), label=c("<0.001", "0.01", "0.1", "1", "10", "100"), cex=0.7)

	mtext(text="time (seconds)", side=2, line=2.5)
	mtext(text="precision", side=1, line=2)
	
	
	txt = paste0("Work-precision diagram, ", k+1, " areas, ", 2^(k+1), " states")
	title(txt)
	
	abline(h=1, lty="dashed", col="grey")
	
	legend(x="topleft", legend=solver_names[-1][keep_cols_TF], col=colors[keep_cols_TF], lty="solid", lwd=1, text.col=colors[keep_cols_TF], pch=(1:(ncol(tmpmat)-1))[keep_cols_TF])
	
	for (j in 2:ncol(tmpmat))
		{
		y = tmpmat[,j]
		
		if (all(is.nan(y)) == TRUE)
			{
			next()
			}
		
		TF = is.nan(y)
		y = y[TF == FALSE]
		x = abstol[TF == FALSE]
		lines(x, y, col=colors[j-1])
		points(x, y, pch=j-1, col=colors[j-1])
		}
	}

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


