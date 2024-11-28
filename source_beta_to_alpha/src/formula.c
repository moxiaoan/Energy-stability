#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include "ScLETD.h"


double CH_f(int n)
{
	double f;
	double ett1, ett2, c2;
	double pfeta1, qfeta1, pfeta2, qfeta2, dGBCC_1, dGBCC_2, dGHCP_1, dGHCP_2;
	
	//ett1 = 1.0 - ac[0].u;
	//ett2 = 1.0 - ac[1].u;
	c2 = 1.0 - ch[0].u - ch[1].u;
	if (ch[0].u < 0.0 || ch[1].u < 0.0)
	{
		//printf("bad!\n");
	}
			
	dGBCC_1=-GBCCTI+GBCCAL+R_a*T_a*log(ch[0].u/c2)+ch[1].u*(BL12_0+BL12_1*(ch[0].u-ch[1].u))+ch[0].u*ch[1].u*BL12_1-ch[1].u*(BL32_0+BL32_1*(c2-ch[1].u)+BL32_2*(c2-ch[1].u)*(c2-ch[1].u))-c2*ch[1].u*(BL32_1+2.0*BL32_2*(c2-ch[1].u))+(c2-ch[0].u)*(BL13_0+BL13_1*(ch[0].u-c2)+BL13_2*(ch[0].u-c2)*(ch[0].u-c2))+ch[0].u*c2*(2.0*BL13_1+4*BL13_2*(ch[0].u-c2))+BL132_0*ch[1].u*(c2-ch[0].u);
           
        dGBCC_2=-GBCCTI+GHSERV+R_a*T_a*log(ch[1].u/c2)+ch[0].u*(BL12_0+BL12_1*(ch[0].u-ch[1].u))-ch[0].u*ch[1].u*BL12_1+(c2-ch[1].u)*(BL32_0+BL32_1*(c2-ch[1].u)+BL32_2*(c2-ch[1].u)*(c2-ch[1].u))+c2*ch[1].u*(-2.0*BL32_1-4*BL32_2*(c2-ch[1].u))-ch[0].u*(BL13_0+BL13_1*(ch[0].u-c2)+BL13_2*(ch[0].u-c2)*(ch[0].u-c2))+ch[0].u*c2*(BL13_1+2.0*BL13_2*(ch[0].u-c2))+BL132_0*ch[0].u*(c2-ch[1].u);

	dGHCP_1=-GHSERTI+GHCPAL+R_a*T_a*log(ch[0].u/c2)+ch[1].u*(HL12_0+HL12_1*(ch[0].u-ch[1].u))+ch[0].u*ch[1].u*HL12_1-ch[1].u*HL32_0+(c2-ch[0].u)*(HL13_0+HL13_1*(ch[0].u-c2)+HL13_2*(ch[0].u-c2)*(ch[0].u-c2))+ch[0].u*c2*(2.0*HL13_1+4*HL13_2*(ch[0].u-c2))+ch[1].u*(c2-ch[0].u)*(ch[0].u*HL132_0+c2*HL132_1+ch[1].u*HL132_2)+ch[0].u*ch[1].u*c2*(HL132_0-HL132_1);
            
        dGHCP_2=-GHSERTI+GHCPV+R_a*T_a*log(ch[1].u/c2)+ch[0].u*(HL12_0+HL12_1*(ch[0].u-ch[1].u))-ch[0].u*ch[1].u*HL12_1+(c2-ch[1].u)*HL32_0-ch[0].u*(HL13_0+HL13_1*(ch[0].u-c2)+HL13_2*(ch[0].u-c2)*(ch[0].u-c2))+ch[0].u*c2*(HL13_1+2.0*HL13_2*(ch[0].u-c2))+ch[0].u*(c2-ch[1].u)*(ch[0].u*HL132_0+c2*HL132_1+ch[1].u*HL132_2)+ch[0].u*ch[1].u*c2*(HL132_2-HL132_1);
	
	//pfeta1=ac[0].u*ac[0].u*ac[0].u*(10.0-15.0*ac[0].u+6.0*ac[0].u*ac[0].u) + ac[1].u*ac[1].u*ac[1].u*(10.0-15.0*ac[1].u+6.0*ac[1].u*ac[1].u);
        //pfeta2=1 - (ac[0].u*ac[0].u*ac[0].u*(10.0-15.0*ac[0].u+6.0*ac[0].u*ac[0].u) + ac[1].u*ac[1].u*ac[1].u*(10.0-15.0*ac[1].u+6.0*ac[1].u*ac[1].u));

	//pfeta2=ac[0].u*ac[0].u*ac[0].u*(10.0-15.0*ac[0].u+6.0*ac[0].u*ac[0].u)+ac[1].u*ac[1].u*ac[1].u*(10.0-15.0*ac[1].u+6.0*ac[1].u*ac[1].u);
        //pfeta1=ett1*ett1*ett1*(10.0-15.0*ett1+6.0*ett1*ett1)+ett2*ett2*ett2*(10.0-15.0*ett2+6.0*ett2*ett2);
	
	pfeta2=ac[0].u*ac[0].u*ac[0].u*(10.0-15.0*ac[0].u+6.0*ac[0].u*ac[0].u)+ac[1].u*ac[1].u*ac[1].u*(10.0-15.0*ac[1].u+6.0*ac[1].u*ac[1].u)+ac[2].u*ac[2].u*ac[2].u*(10.0-15.0*ac[2].u+6.0*ac[2].u*ac[2].u)+ac[3].u*ac[3].u*ac[3].u*(10.0-15.0*ac[3].u+6.0*ac[3].u*ac[3].u);
        pfeta1=1.0-pfeta2;
	switch (n)
	{
		case 0:
		{
			f = (pfeta1*dGBCC_1+pfeta2*dGHCP_1)/gnormal;
			break;
		}
		case 1:
		{
			f = (pfeta1*dGBCC_2+pfeta2*dGHCP_2)/gnormal;
			break;
		}
	}	


	return f;
}
