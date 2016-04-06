#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<cstdlib>

using namespace std;

const int transitiontable[4]={2,3,0,1}; // C<-> T   A<->G
const char AAcode[4]={'A','C','G','T'};
const int complement[4]={3,2,1,0};                // A <-> T   C <-> G
const int circle=1000;
const int phred=33;

struct mitpos{
	char std,snp; //reference base
};

inline int multidx(int i, int j,int k, int l,int readln){  // a combination of the following: {NA,NC,NG,NT,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT }
	k+=(j+1)*4;
	return i*20*readln+k*readln+l;
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


int readref(vector<mitpos> & mitochondrial, ifstream & in){
	mitpos tmp;
	char buf[8192];
	int i,j;
	in.getline(buf,8100);
	tmp.snp=-1;
	while(1){
		in>>tmp.std;
		if(in.eof()) break;
		tmp.std=decode(tmp.std);
		mitochondrial.push_back(tmp);
	}
	return 0;
}


int readvcf(vector<mitpos> & mitochondrial, ifstream & in){
        char buf[81920];
        int i,pos;
        //head skip
        in>>buf>>buf>>buf;
        in.getline(buf,81000);
        if(buf[0]!='\t'){
                cerr<<"No alt allele shown in hsd file\n";
            //    for(i=0;i< mitochondrial.size();i++)  mitochondrial[i].snp=-1;
                return 0;
        }
        pos=0;
        for(i=1;buf[i];i++){
                if(buf[i]>='0' && buf[i]<='9') { pos=pos*10+buf[i]-'0'; continue; }
                if(buf[i]=='A' || buf[i]=='C' || buf[i]=='G' || buf[i]=='T') {
                        mitochondrial[pos-1].snp=decode(buf[i]);
                        if(mitochondrial[pos-1].snp == mitochondrial[pos-1].std) mitochondrial[pos-1].snp=-1;
                }
                for(;buf[i] && buf[i]!='\t' ;i++);
                if(buf[i]==0) break;
                pos=0;
        }
        return 0;
}

/*
int readvcf(vector<mitpos> & mitochondrial, ifstream & in){
	char buf[81920];
	int i,j,pos;
	char ref,alt;
	
	while(1){
		in.getline(buf,81000);
		if(in.eof()) break;
		if(buf[0]=='#') continue;
		for(i=0;buf[i] && buf[i]!='\t';i++);
		for(i++,pos=0;buf[i]>='0' && buf[i]<='9';i++) pos=pos*10+buf[i]-'0';
		for(i++;buf[i] && buf[i]!='\t';i++);
		i++; ref=buf[i]; alt=buf[i+2];
		if(decode(ref)!=mitochondrial[pos-1].std){ 
			cerr<<buf<<endl;
			cerr<<"Error in vcf file : "<<ref<<"\t"<<AAcode[mitochondrial[pos-1].std]<<endl; return -1; 
		}
		mitochondrial[pos-1].snp=decode(alt);
		if(mitochondrial[pos-1].snp==-1) { 
			cerr<<buf<<endl;
			cerr<<"Error in vcf file : "<<alt<<endl; return -1; 
		}
		if(pos<=circle) mitochondrial[mitochondrial.size()-circle+pos-1].snp=mitochondrial[pos-1].snp;
	}
	return 0;
}
*/
int readsam(vector<mitpos> & mitochondrial, ifstream & in, int* err, int* all, int readln ){
	char qname[512],qname2[512],rname[512],rname2[512],cigar[512],cigar2[512],rnext[512],rnext2[512],seq[8192],seq2[8192],qual[8192],qual2[8192],addfield[32768],addfield2[32768];
	int flag,flag2,pos,pos2,mapq,mapq2,pnext,pnext2,tlen,tlen2,i,j,k,creads=0,aa,pa,ma,matchlow,matchhigh;
	
	while(1){
		in>>qname>>flag>>rname>>pos>>mapq>>cigar>>rnext>>pnext>>tlen>>seq>>qual;
		if(in.eof()) break;
		in.getline(addfield,32700);
		if((flag&0xF)!=0x3) continue;
		creads++;
		if(creads%100000==0) cerr<<creads<<" Reading "<<qname<<endl;
	//	if(flag&0x10) reverse(seq,qual);
		if((flag&0xC0)==0x40) k=readln/2; else k=0; 
		i=0; // position in read
		mapq2=0; // current segment in cigar
		pos2=pos-1; // current pos in genome
		for(j=0;cigar[j];j++){ 
			switch(cigar[j]){
				case 'D':
/*						pnext2=0;
						aa=decode(seq[i+pnext2]);
						if(aa==-1)      continue;
						if(flag&0x10){
							aa=complement[aa];
							tlen2=readln/2-1-i-pnext2; //tlen2 cycle, not current pos
							if(tlen2>0) pa=decode(seq[i+pnext2+1]);
							else pa=-1;
							if(pa!=-1) pa=complement[pa];
						}
						else {
							tlen2=i+pnext2;
							if(tlen2>0) pa=decode(seq[i+pnext2-1]);
							else pa=-1;
						}
						
						all[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
						err[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
						if(pa!=-1){ 
							all[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
							err[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
						}
*/						
					
					pos2+=mapq2;
					mapq2=0;
					break;
				case 'M':
					for(matchlow=j-1; matchlow>=0 && cigar[matchlow]>='0' && cigar[matchlow]<='9'; matchlow--);
					if(matchlow==-1) pnext2=0; else pnext2=5;
					if(cigar[j+1]==0) matchhigh=mapq2; else matchhigh=mapq2-5;
					for(;pnext2<matchhigh;pnext2++){
						ma=aa=decode(seq[i+pnext2]);
						if(aa==-1)      continue;
						if(flag&0x10){
							aa=complement[aa];
							tlen2=readln/2-1-i-pnext2; //tlen2 cycle, not current pos
							if(tlen2>0) pa=decode(seq[i+pnext2+1]);
							else pa=-1;
							if(pa!=-1) pa =complement[pa];
						}
						else {
							tlen2=i+pnext2;
							if(tlen2>0) pa=decode(seq[i+pnext2-1]);
							else pa=-1;
						}
					//	if(tlen2<0 || tlen2>=125) cerr<<"Error in tlen2\n";
						all[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
						if(pa!=-1) all[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
						if(ma != mitochondrial[pos2+pnext2].std && ma != mitochondrial[pos2+pnext2].snp){ 
							err[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
							if(pa!=-1) err[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
						}
					}
					pos2+=mapq2;
					i+=mapq2;
					mapq2=0;
					break;
				case 'I':
/*					for(pnext2=0;pnext2<mapq2;pnext2++){
						aa=decode(seq[i+pnext2]);
						if(aa==-1)      continue;
						if(flag&0x10){
							aa=complement[aa];
							tlen2=readln/2-1-i-pnext2; //tlen2 cycle, not current pos
							if(tlen2>0) pa=decode(seq[i+pnext2+1]);
							else pa=-1;
							if(pa!=-1) pa = complement[pa];
						}
						else {
							tlen2=i+pnext2;
							if(tlen2>0) pa=decode(seq[i+pnext2-1]);
							else pa=-1;
						}
						
						all[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
						err[multidx(qual[i+pnext2]-phred,pa,aa,tlen2+k,readln)]++;
						if(pa!=-1){ 
							all[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
							err[multidx(qual[i+pnext2]-phred,-1,aa,tlen2+k,readln)]++;
						}
					} */
					i+=mapq2;
					mapq2=0;
					break;
				case 'S':
				case 'H':	
					i+=mapq2;
					mapq2=0;
					break;
				default: mapq2=mapq2*10+cigar[j]-'0';
			}
		}
	}
			
	
	return 0;
}



int main(int argc, char** argv){
	ifstream ref(argv[1]),samfile(argv[2]),vcffile(argv[3]);
	ofstream calitable(argv[4]);
	vector<mitpos> refgenome;
	int readln=atoi(argv[5]);
	readln=2*readln;
	int i,j,k,l,tiv, *err, *all,n,allbase,errbase,bk,basecount[60],bcount[60],bi,bestbk,n1,n2; 
	double *escore; // qual,aa,direction,pos,
	double bq,mx,my,cov,var,syy,beta1,alpha1,beta2,alpha2,finalb1,finalb2,finala1,finala2,bestse;
//	if(readln<1 || readln >500) { cerr<<"readln does not correct\n"; return -1; }
	cerr<<"Readln = "<<readln<<endl;
	
	err = new int [60*20*readln];
	all = new int [60*20*readln]; // a combination of the following: {NA,NC,NG,NT,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT }
	escore = new double [60*20*readln]; //phred score is set up to 59
	

	for(i=0;i<60;i++)
                
           for(k=0;k<20;k++)		
				for(l=0;l<readln;l++){ 
//					cerr<<i<<"\t"<<j<<"\t"<<k<<"\t"<<l<<endl;
					escore[multidx(i,k/4-1,k%4,l,readln)]=0;
					if(i==2){
						err[multidx(i,k/4-1,k%4,l,readln)]=6310;
                                        	all[multidx(i,k/4-1,k%4,l,readln)]=10000;
					}
					else{
						err[multidx(i,k/4-1,k%4,l,readln)]=0;
                                        	all[multidx(i,k/4-1,k%4,l,readln)]=0;
					}
				}					
			
		
	
	cerr<<"Done with allocation\n";

	readref(refgenome,ref);
	cerr<<"Done with reference\n";
	ref.close();
	readvcf(refgenome,vcffile);
	cerr<<"Done with vcf\n";
	vcffile.close();
	readsam(refgenome,samfile,err,all,readln); //in paired mode
	cerr<<"Done with reads\n";
	cerr<<"Ref size = "<<refgenome.size()<<endl;
	
	
		for(k=0;k<20;k++)
			for(l=0;l<readln;l++){
				n=0;  bi=0;
				for(i=0;i<60;i++){
				
					allbase=all[multidx(i,k/4-1,k%4,l,readln)]/100; errbase=err[multidx(i,k/4-1,k%4,l,readln)];
					//if(allbase>0) basecount[i]=1; else basecount[i]=0;
					
					if(allbase==0 || errbase==0){ basecount[i]=0;  continue; }
					basecount[i]=allbase;
					bcount[bi]=i; bi++;
					n+=allbase;
			
					
					bq=(-10)*log10((double)(errbase)/(all[multidx(i,k/4-1,k%4,l,readln)]));
					
					escore[multidx(i,k/4-1,k%4,l,readln)]=bq;
						
				}
				bestbk=-1; finala1=finalb1=0;
					
				 finala2=finalb2=0; bestse=1e300; 

				if(bi>5){ // need at least 6 samples to have the segment regression
					for(bk=3;bk<bi-2;bk++){
						n1=n2=0; 
						mx=my=cov=var=syy=0;
						if(bcount[bk]>25) {
							bk=1; beta1=alpha1=0;
							for(i=bcount[bk];i<=bcount[bi-1];i++){
								n2+=basecount[i];
								mx+=i*basecount[i];
								my+=escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
								cov+=i*escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
								var+=i*i*basecount[i];
							}
							beta2=(n2*cov-mx*my)/(n2*var-mx*mx);
							alpha2=my/n2-beta2*mx/n2;
							for(i=bcount[bk];i<=bcount[bi-1];i++){
								bq=escore[multidx(i,k/4-1,k%4,l,readln)]-beta2*i-alpha2;
								syy+=basecount[i]*bq*bq;							
							}
							if(syy<bestse){ bestse=syy; bestbk=bcount[bk]; finala1=alpha1; finala2=alpha2; finalb1=beta1; finalb2=beta2;}
							break;
						}
						for(i=0;i<bcount[bk];i++){
							n1+=basecount[i];
							mx+=i*basecount[i];
							my+=escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
							cov+=i*escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
							var+=i*i*basecount[i];
						}
						beta1=(n1*cov-mx*my)/(n1*var-mx*mx);
						alpha1=my/n1-beta1*mx/n1;
						for(i=0;i<bcount[bk];i++){
							bq=escore[multidx(i,k/4-1,k%4,l,readln)]-beta1*i-alpha1;
							syy+=basecount[i]*bq*bq;							
						}
						//the first segment is done
						mx=my=cov=var=0;
					//	if(bcount[bk]>17) i=17; else i=bcount[bk];
						for(i=bcount[bk];i<=bcount[bi-1];i++){
							n2+=basecount[i];
							mx+=i*basecount[i];
							my+=escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
							cov+=i*escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
							var+=i*i*basecount[i];
						}
						beta2=(n2*cov-mx*my)/(n2*var-mx*mx);
						alpha2=my/n2-beta2*mx/n2;
						if(bcount[bk]>20) i=20; else i=bcount[bk];
						for(;i<=bcount[bi-1];i++){
							bq=escore[multidx(i,k/4-1,k%4,l,readln)]-beta2*i-alpha2;
							syy+=basecount[i]*bq*bq;							
						}
						if(syy<bestse){ bestse=syy; bestbk=bcount[bk]; finala1=alpha1; finala2=alpha2; finalb1=beta1; finalb2=beta2;
				//				for(tiv=bestbk-1;tiv>bcount[bk-1] && finalb2*tiv+finala2<finalb1*tiv+finala1;tiv--);
				//				bestbk=tiv+1;
						}
					}
									
				}
				else if(bi>2){
					n2=0; bk=1; beta1=alpha1=0;
					mx=my=cov=var=syy=0;
					for(i=bcount[bk];i<=bcount[bi-1];i++){
						n2+=basecount[i];
						mx+=i*basecount[i];
						my+=escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
						cov+=i*escore[multidx(i,k/4-1,k%4,l,readln)]*basecount[i];
						var+=i*i*basecount[i];
					}
						beta2=(n2*cov-mx*my)/(n2*var-mx*mx);
						alpha2=my/n2-beta2*mx/n2;
						bestbk=bcount[bk]; finala1=alpha1; finala2=alpha2; finalb1=beta1; finalb2=beta2;
				//		for(tiv=bestbk-1;tiv>bcount[bk-1] && finalb2*tiv+finala2<finalb1*tiv+finala1;tiv--);
				//				bestbk=tiv+1;
				}
				calitable<<bestbk<<"\t"<<finala1<<"\t"<<finalb1<<"\t"<<finala2<<"\t"<<finalb2<<"\t";
//				if(k/4==0) calitable<<'N'; else calitable<<AAcode[k/4-1];
//				calitable<<AAcode[k%4]<<"\tCycle ";
//				if(l>=125) calitable<<"-"<<l-125+1<<endl;
//				else calitable<<"+"<<l+1<<endl;
				for(i=0;i<60;i++) calitable<<escore[multidx(i,k/4-1,k%4,l,readln)]<<"\t";
				calitable<<endl;
//				for(i=0;i<60;i++) calitable<<err[multidx(i,k/4-1,k%4,l,readln)]<<"\t";
//				calitable<<endl;
//				for(i=0;i<60;i++) calitable<<all[multidx(i,k/4-1,k%4,l,readln)]<<"\t";
//				calitable<<endl;
			}
	
	

	
	return 0;
}
