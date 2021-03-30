function weight=parker_weight(beta,alpha_cor,alpha_thre)
[cor_beta,cor_alpha]=meshgrid(beta,alpha_cor);
cor_beta=cor_beta';cor_alpha=cor_alpha';
alpha_1=2*(alpha_thre+cor_alpha);
alpha_2=pi+2*cor_alpha;
alpha_3=pi+2*alpha_thre;
id1=(cor_beta>=0)&(cor_beta<alpha_1);
id2=(cor_beta>=alpha_1)&(cor_beta<alpha_2);
id3=(cor_beta>=alpha_2)&(cor_beta<alpha_3);
id4=(cor_beta>=alpha_3)&(cor_beta<=2*pi);
assert(sum(id1+id2+id3+id4==1,'all')==size(id1,1)*size(id1,2));
weight=zeros(size(id1));
tmp1=(sin(pi/4*cor_beta./(alpha_thre+cor_alpha))).^2;
weight(id1)=tmp1(id1);
weight(id2)=1;
tmp3=(sin(pi/4*(pi+2*alpha_thre-cor_beta)./(alpha_thre-cor_alpha))).^2;  %%%dissertion has a error here
weight(id3)=tmp3(id3);
weight(id4)=0;
end