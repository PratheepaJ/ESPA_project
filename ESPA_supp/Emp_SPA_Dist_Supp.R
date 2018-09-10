#Source the function to compute the saddlepoint pdf
Emp_SPA_Dist_Supp=function(time, status){
      #time=observed time vector
      
      #status=indicate whether observation is censored (0) or uncensored (1)
      
      #To avoid negative values for a log transformation, add 1 to each observed time
      
      #log transformation
      Y=log(time+1)
      
      # sort the times, and the status keeping it with the corresp time (concomitant) 
      delta=status[sort.list(Y)]; z=sort(Y); n=length(z);
      
      #############################################################
      # Compute fhat(Z(i)) when delta(i)=1 :call it as pw[i]      #
      #############################################################
      pw=rep(0,n); pw[1]=delta[1]/n; prod.term=1;
      for (i in 2:n){
            prod.term = prod.term * (1 - delta[i-1]/(n-i+2))
            pw[i] = prod.term * delta[i]/(n-i+1)
      } 
            
      #############################################################
      # Compute phi_n                                             #
      #############################################################
      p.str=sum(pw);
      
      # decide whether or not to exp tail completion based on delta[n]
      phi=0; if (delta[n]==0) {phi = -log(1-p.str)/z[n]};
      
      #############################################################
      # Empirical CGF of Z                                        #
      #############################################################
      
      K.cgf=function(s){
            Mds=sum(pw*exp(s*z)); Mcs=0;
            if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
            log(Mds+Mcs)
      }
      
      #############################################################
      # Second derivative of Empirical CGF of Z                   #
      #############################################################
      
      K2.exp=function(s){
            Mds=sum(pw*exp(s*z)); Mcs=0;
            if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
            Ms=Mds+Mcs;
            #   first primes
            Mds.prime=sum(pw*z*exp(s*z)); Mcs.prime=0;
            if (delta[n]==0) {Mcs.prime=Mcs*(z[n]+1/(phi-s))};
            Ms.prime=Mds.prime+Mcs.prime
            #   double primes
            Mds.prime2=sum(pw*z^2*exp(s*z)); Mcs.prime2=0;
            if (delta[n]==0) {Mcs.prime2=Mcs.prime*(z[n]+1/(phi-s))+Mcs/(phi-s)^2};
            Ms.prime2=Mds.prime2+Mcs.prime2
            #   put it all together
            K.prime2=(Ms*Ms.prime2-Ms.prime^2)/Ms^2
            K.prime2
      }
      
      #############################################################
      # Adjusted CGF to optimze (minimize): to find the sh        #
      #############################################################
      
      K.adj.exp=function(s,t){
            Mds=sum(pw*exp(s*z)); Mcs=0;
            
            if (delta[n]==0){
                  Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
            }
            Ms=Mds+Mcs
            
            if(round(Ms,3)==0|round(K2.exp(s),3)<=0){
                  return(10000)
            }
            else{
                  return(K.cgf(s)-t*s)
            }
      }
      
      #############################################################
      #       Use minimization to find the saddlepoint            #
      #############################################################
      saddle.me<-function(t){
            suppressWarnings(optim(0,K.adj.exp,t=t,method="Nelder-Mead",hessian=FALSE)$par)
      }
      
      ############################################################
      #  find the lower support of the failure time              #
      #  saddelpoint should be less than or equal to phi         #
      ############################################################
      
      M.mgf=function(s){
            Mds=sum(pw*exp(s*z)); Mcs=0;
            if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
            Ms=Mds+Mcs;
            return(Ms)
      }
      
      x.min=min(z)
      sh.mi=saddle.me(x.min)
      
      while(sh.mi<phi & (exp(x.min)-1)>0){#for small values of time, s is negative
            x.min=x.min-.01
            sh.mi=saddle.me(x.min)
            if(round(M.mgf(sh.mi),3)<=0|round(K2.exp(sh.mi),3)<=0){break}
      }  
      
      ############################################################
      #  Determine the max of the support:                       #
      ############################################################
      
      x.max=max(z)
      pdf=Emp_SPA_pdf_t(exp(x.max)-1,time,status)
      sh.m=saddle.me(x.max)
      if(delta[n]==0){
            while(sh.m<phi & round(pdf,3)>0){
                  x.max=x.max+.01#it will increase by .01 but the logic in the while fail thereafter
                  pdf=Emp_SPA_pdf_t(exp(x.max)-1,time,status)
                  sh.m=saddle.me(x.max)
            }  
      }else{
            while(round(pdf,3)>0){
                  x.max=x.max+.01#it will increase by .01 but the logic in the while fail thereafter
                  pdf=Emp_SPA_pdf_t(exp(x.max)-1,time,status)
                  sh.m=saddle.me(x.max)
            } 
      }
      
      
      minSp=round(exp(x.min+.01)-1,3)
      
      maxSp=round(exp(x.max-.01)-1,3)
      
      return(c(minSp,maxSp))
      
}
