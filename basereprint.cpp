#include<iostream>
#include<fstream>
#include<cstdlib>

using namespace std;

const int phred = 33;
const char AAcode[4]={'A','C','G','T'};
const int complement[4]={3,2,1,0};  

struct regression{
	int bk;
	double a1,b1,a2,b2;
	int cap;
};

inline int multidx(int j,int k, int l,int readln){  // a combination of the following: {NA,NC,NG,NT,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT }
	k+=(j+1)*4;
	return k*readln+l;
}

int decode(char ch){
	int val;
	switch(ch){
				case 'A':
				case 'a': val=0; break;
				case 'C':
				case 'c': val=1; break;
				case 'G':
				case 'g': val=2; break;
				case 'T':
				case 't': val=3; break;
				default: val=-1;
	}
	return val;
}

int main(int argc, char** argv){
	char qname[512],rname[512],cigar[512],rnext[512],seq[8192],qual[8192],addfield[32768];
	int flag,pos,mapq,pnext,tlen,readln,aa,pa,cycle,rev,k,empq;
	double sampleph[60];
	int capdecision[60];
	regression *table;
	int i,j,l,overallcap=0;
	ifstream recali(argv[1]);
	readln=2*atoi(argv[2]);
	table= new regression[readln*20];
	for(i=0;i<readln*20;i++){
		recali>>table[i].bk>>table[i].a1>>table[i].b1>>table[i].a2>>table[i].b2;
		for(k=j=0;j<60;j++){
			recali>>sampleph[k];
		//	cerr<<capdecision[0][j]<<"\t";
			if(sampleph[k]>2){
				capdecision[k]=j;
				k++;
				if(j>overallcap) overallcap=j;
			}
		}
		recali.getline(addfield,30000);
		for(k--;k>0;k--) if(sampleph[k]>sampleph[k-1]) break;
		table[i].cap=capdecision[k];
	//	cerr<<endl;
	}
		
	while(1){
		cin>>addfield[0];
		if(addfield[0]!='@'){ cin.putback(addfield[0]); break; }
		cin.getline(addfield+1,30000);
		cout<<addfield<<endl;
	}
	
	// header is copied;
	while(1){
		cin>>qname>>flag>>rname>>pos>>mapq>>cigar>>rnext>>pnext>>tlen>>seq>>qual;
		if(cin.eof()) break;		
		cin.getline(addfield,32768);
		if((flag&0x1)==0x0) continue;
	//	if(flag&0x1){ // worked only on original reads
			if((flag&0xC0)==0x40) rev=readln/2; else rev=0; 
			for(i=0;i<readln/2;i++){
				if(qual[i]<phred+5) continue; //low quality in raw qual is ignored
				aa=decode(seq[i]);
				if(flag&0x10){
					aa=complement[aa];
					cycle=readln/2-1-i; 
					if(cycle>0) pa=decode(seq[i+1]);
					else pa=-1;
					if(pa!=-1) pa =complement[pa];
				}
				else {
					cycle=i;
					if(cycle>0) pa=decode(seq[i-1]);
					else pa=-1;
				}
				k=multidx(pa,aa,cycle+rev,readln);
				if(table[k].bk!=-1 && pa!=-1) k=multidx(-1,aa,cycle+rev,readln);
				if(table[k].bk!=-1){
					if(qual[i]-phred<table[k].bk) empq=table[k].b1*(qual[i]-phred)+table[k].a1;
					else  empq=table[k].b2*(qual[i]-phred)+table[k].a2;
					l=table[k].b2*(table[k].cap)+table[k].a2;
					if(l>overallcap) l=overallcap;
					if(empq<2) empq=2;
					if(empq>l) empq=l;
					qual[i]=phred+empq;
				}
			}
	//	}
		cout<<qname<<'\t'<<flag<<'\t'<<rname<<'\t'<<pos<<'\t'<<mapq<<'\t'<<cigar<<'\t'<<rnext<<'\t'<<pnext<<'\t'<<tlen<<'\t'<<seq<<'\t'<<qual<<addfield<<endl;
	}
}
