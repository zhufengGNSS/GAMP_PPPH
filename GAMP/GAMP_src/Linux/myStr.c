#include "gamp.h"

//convert number to string
extern void num2str(int num, char *str, int len)
{
	int i,j=0,n;
	char tmp[MAXCHARS]={'\0'};

	n=len-(int)(log10((float)num))-1;
	if (n<0) {
		tmp[0]='\0';
		return;
	}

	for (i=0;i<n;i++) {
		strcpy(&tmp[j++],"0");
	}
	sprintf(str,"%s%d%c",tmp,num,'\0');
}

//string clipping
extern void xStrMid (char *szDest, const int nPos, const int nCount, char *szSrc)
{
	int i,n;
	char *str,c;

	n=strlen(szSrc);
	if (n<=0) return;

	str=szSrc+nPos;
	for (i=0;i<nCount;i++) {
		c=*(str+i);
		if (c) {
			*(szDest+i)=c;
		}
		else {
			*(szDest+i)='\0';
			break;
		}
	}
	*(szDest+nCount)	='\0';
}

extern void trimSpace(char *strsrc)
{
    int i=0,j=0,ps,pe;
    int len=strlen(strsrc);
    char str[MAXCHARS+1];
    
    if (len<=0) return;
    
    str[0]='\0';
    strcpy(str,strsrc);
    
    ps=0;
    for (i=0;i<len;i++) {
        if (*(str+i)!=' '&&*(str+i)!='\t') {
            ps=i;
            break;
        }
    }

    pe=ps;
    for (j=len-1;j>=0;j--) {
        if (*(str+j)!=' '&&*(str+j)!='\t'&&*(str+j)!='\n') {
            pe=j;
            break;
        }
    }

    if (pe==ps) 
        *(str+pe)='\0';
    else
        *(str+pe+1)='\0';
    
    strcpy(strsrc,str+ps);
}

extern void cutFilePathSep(char *strPath)
{
	int i,len;

	for (i=0;i<4;i++) {
		len=strlen(strPath);

		if (len<=0) 
			break;

		if (strPath[len-1]==FILEPATHSEP) 
			strPath[len-1]='\0';
		else
			break;
	}
}
