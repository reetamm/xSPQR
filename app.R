library(shiny)
library(shinydashboard)
create_plot = function(K,seed,xi,p.a,p.b,c1,c2){
  set.seed(seed)
  pred <- LaplacesDemon::rdirichlet(1,alpha=rep(1,K))
  
  library(splines2)
  library(ggplot2)
  basis <- function(y,K,integral=FALSE){
    library(splines2)
    knots <- seq(1/(K-2), 1-1/(K-2), length=K-3)
    B     <- mSpline(y, knots = knots, 
                     Boundary.knots=c(0,1),
                     intercept=TRUE, degree = 2,
                     integral=integral)
    return(t(B))}
  
  spqrk <- function(pred,Y=NULL,nY=101,tau=0.5,type='QF'){
    # pred        <- as.matrix(model(X))
    n <- nrow(pred)
    ntau = length(tau)
    n.knots <- ncol(pred)
    if(is.null(Y) | type=='QF')
      Y <- seq(0,1,length.out = nY)
    B <- (basis(Y , n.knots,integral = (type!='PDF')))
    if(ncol(B)!=n)
      df1  <- pred%*%B
    if(ncol(B)==n)
      df1  <- colSums(B*t(pred))
    if(type!='QF'){
      return(df1)
    }
    if(type=='QF'){
      qf1 = matrix(NA,n,ntau)
      for(i in 1:n)
        qf1[i,] <- stats::approx(df1[i,], Y, xout=tau, ties = list("ordered", min))$y
      return(qf1)
    }
  }


  F.SPQR=function(y){
    y[y<0]=0
    y[y>1]=1
    spqrk(pred,Y=y,type = 'CDF')  
  }
  f.SPQR=function(y){
    try(out<-spqrk(pred,Y=y,type = 'PDF') , silent=T)
    out[y<0]=0
    out[y>1]=0
    out
  }
  Finv.SPQR=function(y){
    spqrk(pred,tau =y,type = 'QF')  
    
    
  }
  F.GPD=function(y,u=0,scale=1,shape=0.1){
    evd::pgpd(y,u,scale,shape)
  }
  f.GPD=function(y,u=0,scale=1,shape=0.1){
    try(out<-(1/scale)*(1+shape*(y-u)/scale)^(-1/shape-1), silent=T)
    out[1+shape*(y-u)/scale<=0]=0
    out
  }
  
  F.blend=function(y,F.S,F.G,p){
    
    F.S^(1-p)*F.G^(p)
  }
  weight=function(y,a,b,c1=10,c2=10){
    pbeta((y-a)/(b-a),c1,c2)
  }
  weight.prime=function(y,a,b,c1=10,c2=10){
    dbeta((y-a)/(b-a),c1,c2)/(b-a)
  }
  
  
  
  a=Finv.SPQR(p.a)
  a
  b=Finv.SPQR(p.b)
  b
  
  sigma.val=function(a,b,p.a,p.b,xi){
    xi*(a-b)/(((1-p.a)^(-xi))-((1-p.b)^(-xi)))
  }
  u.val=function(a,p.a,sigma,xi){
    a-sigma/xi*((1-p.a)^(-xi)-1)
  }
  
  
  sig<-c(sigma.val(a,b,p.a,p.b,xi))
  
  u<-c(u.val(a,p.a,sig,xi))
  
  y=seq(0,1.5,length=500)
  
  probs=apply(as.matrix(y),1,function(x){
    
    F.blend(x,F.SPQR(x),F.GPD(x,u=u,scale=sig,shape=xi),weight(x,a,b,c1,c2))
  })
  par(mfrow=c(2,1))
  
  
  f.blend=function(y,F.B,F.S,f.S,F.G,f.G,p,p.prime){
    
    if(F.G==0) return(f.S)
    F.B*(p.prime*log(F.G)+p*f.G/F.G-p.prime*log(F.S)+(1-p)*f.S/F.S)
  }
  
  # F.S=F.SPQR(0.1)
  # F.G=F.GPD(0.1,u=u,scale=sig,shape=xi)
  # p=weight(0.1,a,b)
  # F.B=F.blend(0.1,F.S,F.G,p)
  # 
  # f.S=f.SPQR(0.1)
  # f.G=f.GPD(0.1,u=u,scale=sig,shape=xi)
  # p.prime=weight.prime(0.1,a,b)
  # f.blend(0.1,F.B,F.S,f.S,F.G,f.G,p,p.prime)
  
  
  
  y=seq(0,1.5,length=500)
  # plot(y,f.GPD(y,u,sig,xi),type="l",xlim=c(0,1.5),lwd=3,col="red")
  # points(y,f.SPQR(y),type="l",col="green",lwd=3)
  # abline(v=a,col="orange")
  # abline(v=b,col="orange")
  #points(y,weight.prime(y,a,b),type="l",col="blue")
  
  dens=apply(as.matrix(y),1,function(x){
    F.S=F.SPQR(x)
    F.G=F.GPD(x,u=u,scale=sig,shape=xi)
    p=weight(x,a,b,c1,c2)
    F.B=F.blend(x,F.S,F.G,p)
    
    f.S=f.SPQR(x)
    f.G=f.GPD(x,u=u,scale=sig,shape=xi)
    p.prime=weight.prime(x,a,b,c1, c2)
    if(F.G==0) return(f.S)
    if(F.G>0) return(f.blend(x,F.B,F.S,f.S,F.G,f.G,p,p.prime))
  })
  #points(y,dens,type="l",col="black",lwd=2)
  print(min(dens,na.rm=T))
  #par(mfrow=c(1,1))
  plot(y,F.GPD(y,u,sig,xi),type="l",ylim=c(0,1), xlim=c(0,1.1),
       lwd=3,col="red", ylab = expression(H(paste(y,"|", W), xi)),
       main ="")
  points(y,F.SPQR(y),type="l",ylim=c(0,1),col="green",lwd=3)
  abline(v=a,col="orange")
  abline(v=b,col="orange")
  points(y,probs,type="l",ylim=c(0,1),col="black",lwd=2)
  
  
  plot(y,f.GPD(y,u,sig,xi),type="l",xlim=c(0,1.1),
       ylim=range(dens,na.rm=T),lwd=3,col="red",
       ylab = expression(h(paste(y,"|", W), xi)),
       main = "")
  points(y,f.SPQR(y),type="l",col="green",lwd=3)
  abline(v=a,col="orange")
  abline(v=b,col="orange")
  points(y,dens,type="l",ylim=c(0,1),col="black",lwd=2)
}

# create_plot(K=10,seed = 1, xi, p.a, p.b, c1=20, c2) 


 
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar( sliderInput("sliderxi",withMathJax("$$\\xi$$"), 
                                min=-1, max=1, step=0.025, value=0.1),
                    sliderInput("sliderK",withMathJax("$$K$$"), 
                                min=5, max=50, step=1, value=5),
                    sliderInput("sliderseed",withMathJax("seed"), 
                                min=1, max=5, step=1, value=1),
                    sliderInput("sliderp_a",withMathJax("$$p_a$$"), 
                                min=0.05, max=0.75, step=0.05, value=0.25),
                    sliderInput("sliderp_b",withMathJax("$$p_b$$"), 
                                min=0.75, max=0.95, step=0.05, value=0.8),
                    sliderInput("sliderc_1",withMathJax("$$c_1$$"), 
                                min=5, max=100, step=5, value=25),
                    sliderInput("sliderc_2",withMathJax("$$c_2$$"), 
                                min=5, max=100, step=5, value=5)
                    ),
  dashboardBody(
    fluidRow(column(8,plotOutput('plot')))
  ))

server <- function(input, output, session) { 
  output$plot <- renderPlot({
    create_plot(K=input$sliderK,seed = input$sliderseed, 
                xi = input$sliderxi, p.a = input$sliderp_a, p.b=input$sliderp_b, 
                c1=input$sliderc_1, c2 = input$sliderc_2) 
})
}
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))

