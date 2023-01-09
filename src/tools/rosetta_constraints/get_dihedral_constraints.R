#!/usr/bin/Rscript
require(RSQLite)
require(diffusionMap)
require(rgl)
require(circular)

drv <- dbDriver("SQLite")

dont_plot = TRUE

args<-commandArgs(TRUE)
db_file <- args[1]
cluster_list <- args[2]
outdir <- args[3]


con <- dbConnect(drv, dbname = db_file)
#clusters = read.table(cluster_list)

sql <- "SELECT DISTINCT fullcluster FROM cdr_data WHERE datatag != 'loopKeyNotInPaper'"
clusters = dbGetQuery(con, sql)

plot_dir = paste(sep="/", outdir, "PLOTS")
sink(paste(sep="/", outdir, "MEAN_SD.txt"))
cat("cluster|type|length|n|position|mean_phi|sd_phi|mean_psi|sd_psi|mean_omega|sd_omega\n")
#Loop over each cluster. 
to_circular <- function(data){
  d = circular(data, type="angle", units="degrees")
}
to_degrees <- function(rad){
  d = rad * (180/pi)
}

for (i in clusters$fullcluster){

	sql <- "SELECT length_type,length,dihedrals FROM cdr_data WHERE fullcluster=?"
	d = dbGetQuery(con, sql, toString(i))
	dir = plot_dir
	f = paste(sep="",dir,"/", toString(i), "_KDE.pdf")
	#cat(f)
	#cat("\n")
	if (dont_plot==FALSE){
		pdf(file=f, width=7*5/2, height=7*4/2)
		par(mfrow=c(4,4))
	}
	l = d$length[1]
	#Here, we have to create the dihedral table.
	phi_num = 1; psi_num=2; omega_num = 3
	for (i_length in 1:l){
		#cat(toString(i_length))
		#cat("\n")
		m = matrix(ncol=3, nrow=0)
		for (seq in d$dihedrals){
			#Fix for some reason I didn't grep this out.
			seq = gsub('\r', '', seq)
			seq = gsub('\n', '', seq)
			#cat(seq)
			dihedrals = strsplit(seq, ":")
			#cat(dihedrals[[1]][phi_num])
			#cat("\n")
			phi = as.numeric(dihedrals[[1]][phi_num]); psi = as.numeric(dihedrals[[1]][psi_num]); omega =  as.numeric(dihedrals[[1]][omega_num])
			
			#Get all numbers positive.
			if (phi<0){phi=360+phi}; if (psi<0){psi=360+psi}; if (omega<0){omega=360+omega}
			m <-rbind(m, c(phi,psi,omega))
		}
  	mean_phi = mean(to_circular(m[,1])); sd_phi = to_degrees(sd(to_circular(m[,1])))
		mean_psi = mean(to_circular(m[,2])); sd_psi = to_degrees(sd(to_circular(m[,2])))
		mean_omega = mean(to_circular(m[,3])); sd_omega = to_degrees(sd(to_circular(m[,3])))

		out = paste(sep=" ", toString(i), d$length_type[1], toString(l), toString(length(m[,1])),toString(i_length), toString(mean_phi), toString(sd_phi), toString(mean_psi), toString(sd_psi), toString(mean_omega), toString(sd_omega))
		cat(out)
		cat("\n")

		#To NOT redo plots.
		if (dont_plot){
			phi_num=phi_num+3; psi_num = psi_num+3; omega_num=omega_num+3
			next
		}


		set.seed(1234)
		base = paste(sep="", "Position ", toString(i_length))
		if (length(m[,1])<3){next}
		#PHI
		x <- circular(m[,1], type="angles", units="degrees"); bandwidth <- bw.cv.mse.circular(x)
		head = paste(sep="", "Phi: ", base)
		y=density(x, adjust=1, bw=bandwidth, kernel="vonmises", from=(0), to=2*pi)
		plot(y, shrink=1.75, main=head)
		#PSI
		x <- circular(m[,2], type="angles", units="degrees"); bandwidth <- bw.cv.mse.circular(x)
		head = paste(sep="", "Psi: ", base)
		y=density(x, adjust=1, bw=bandwidth, kernel="vonmises", from=0, to=2*pi)
		plot(y, shrink=1.75, main=head)
		#OMEGA
		x <- circular(m[,3], type="angles", units="degrees"); bandwidth <- bw.cv.mse.circular(x)
		head = paste(sep="", "Omega: ", base)
		y=density(x, adjust=1, bw=bandwidth, kernel="vonmises", from=0, to=2*pi)
		plot(y, shrink=1.75, main=head)

		phi_num=phi_num+3; psi_num = psi_num+3; omega_num=omega_num+3	
			
			
	}

	if (dont_plot==FALSE){
		dev.off()
	}
				
}
#warnings()
sink()	

