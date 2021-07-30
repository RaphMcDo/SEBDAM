#testing SEBDAM
library(devtools)
load_all()
load("test_data_SPA3.RData")
# test_data$Year<-test_data$Year+1
test_data<-subset(test_data,Year>0)
test_data$I<-(test_data$I/1000)/(0.005334*0.8)
test_data$IR<-(test_data$IR/1000)/(0.005334*0.8)
test_growths<-test_growths[-1,]

greg<-setup_mesh(data=test_data,seed=21,model_bound = mod_bound)
greg2<-setup_pred_grid(knots=greg$knots,model_bound=greg$utm_bound)
# test_catch[,23]<-rep(0,length(test_catch[,1]))
# test_catch<-test_catch[,c(23,1:22)]
greg3<-data_setup(data=test_data,growths = test_growths,catch=test_catch,
                  model="SEBDAM",mesh=greg$mesh,bound=greg$utm_bound,obs_mort=F,
                  prior=T,mult_qI = T,spat_approach="spde_aniso",
                  knot_obj=greg$knots,knot_area=greg2$area,separate_R_aniso=F)

greg_fit<-fit_model(greg3)

#If did not work, force aniso R to be aniso B

#Testing TLM
library(devtools)
load_all()
load("test_data_SPA3.RData")
test_data$Year<-test_data$Year+1
test_data$I<-((test_data$I/1000000)/(0.005334*0.8))*3518.7337
test_data$IR<-((test_data$IR/1000000)/(0.005334*0.8))*3518.7337

test_catch[,23]<-rep(0,length(test_catch[,1]))
test_catch<-test_catch[,c(23,1:22)]
tlm_catch<-colSums(test_catch)/1000
greg3<-data_setup(data=test_data,growths=test_growths,catch=tlm_catch,
                  model="TLM",obs_mort=T,prior=T)


greg_fit<-fit_model(greg3)

#Testing Simulations


