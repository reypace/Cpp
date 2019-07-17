/*This is a stat program which takes user input of sample statistics
*and performs all the math necessary to give exact p values for one-sample
*t and z tests, including basic numeric integration. 

Created by Jesse R. Pace
*/

#include <iostream>
#include <cmath>

float integrate(double start, double end); //function to do the integration over the normal distribution
float integrate(double start, double end, double samp); //function to do the integration over the t distribution

int main()
{
	
	float mew;
	float sigmaSq;
	float xbar;
	float lower; 
	float upper;
	float end;  //this will be the end for lower, and the begining for upper
	float pval;	
	float sampsize;
	int tOrZ;	

	std::cout<<"enter sample size"<<"\n";
	std::cin>>sampsize;
	std::cout<<"enter sample mean"<<"\n";
	std::cin>>xbar;
	std::cout<<"enter sample variance"<<"\n";
	std::cin>>sigmaSq;
	std::cout<<"is variance estimated or exact? 1=Est., 2=Exact"<<"\n";
	std::cin>>tOrZ;
	while(tOrZ!=1&tOrZ!=2){
		std::cout<<"please select 1 for Estimated, or 2 for Exact";
		std::cin>>tOrZ;
	}
	std::cout<<"enter mean for null"<<"\n";
	std::cin>>mew;

	/*the probability that xbar is as far or further than the observed xbar away from the null mean, given that the null is true
	* is simply p(XBAR>=xbar)
	*where XBAR is the random var and xbar is the observed (realized value)
	*we simply subtract the mean from each side and divide by std dev to put into std normal form
	*p(XBAR<=|xbar|)=int(-4,xbar-mew/sqrt(sigmaSq))
	*/
	
	if(tOrZ==2){
		end=(xbar-mew)/sqrt(sigmaSq/sampsize);
		end=sqrt(end*end); // first we make sure that end is positive
		upper=integrate(end,4.0);    // plus p xbar is above pos end
		end=-end; //now we look below the neg end
		lower=integrate(-4.0,end);  
	}//end of if ==2
	
	
	
	if(tOrZ==1){
		end=(xbar-mew)/sqrt(sigmaSq/sampsize);
		end=sqrt(end*end); // first we make sure that end is positive
		upper=integrate(end,5.0,sampsize);    // plus p xbar is above pos end
		end=-end; //now we look below the neg end
		lower=integrate(-5.0,end,sampsize);  
	}//end of if ==1
	
	
	pval=lower+upper;
	std::cout<<"\nyour pvalue is";
	std::cout<<pval;
	
	
	return 0;
} // end of int main()



//function to do the integration over the normal distribution
float integrate(double start, double end)
{
	int i;
	double rie;
	double delta;
	double x;
	double expo;
	
	rie=0.0; //manually setting rie to 0, I think C++ might continue to store variables between calls of a function...
	delta=(end-start)/400000;
	i = 1;
	x=start;
	
	while (i<=400000) {
		expo=pow(2.7182818,  ( - (x*x)/2*(1*1) )  ) ;  //c++ we need a seperate function to do powers
		rie= 1/sqrt(2*3.1415926535897*1)*expo*delta+rie;
		i=i+1;
		x=x+delta;
			
	}//end while
	
	return rie;
}//end integrate function

// integration over Student's t
float integrate(double start, double end, double samp)
{
	int i;
	double rie;
	double delta;
	double x;
	double expo;
	double r;
	double constant;
	
	delta=(end-start)/500000;
	rie=0.0;
	r=samp-1;
	i = 1;
	x=start;
	constant = tgamma((r+1)/2)/ (  sqrt(3.1415926535897*r)*tgamma(r/2.0));  //these components are part of the PDF, but there's no sense in looping them in
	
	while (i<=500000) {
		expo=pow( (1+((x*x)/r)), (r+1)/2)   ;  //c++ we need a seperate function to do powers
		rie = (constant/expo)*delta + rie;
		i=i+1;
		x=x+delta;
					
	}//end while
	
	return rie;
}//end integrate function
