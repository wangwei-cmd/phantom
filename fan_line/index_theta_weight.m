function [index]=index_theta_weight(beta_chord,beta)
index=ones(size(beta_chord));
beta_b=beta(1);
assert(beta_b==0);
beta_t=beta(end);
beta_b_chord=beta_chord(1,:,:);
beta_t_chord=beta_chord(end,:,:);
delt_beta=beta(2)-beta(1);
beta_expand=repmat(beta',[1,size(index,2),size(index,3),size(index,4)]);

part1=(beta_chord>=beta_t);
part2=1-part1;
tmp1=beta_expand-beta_b;
tmp2=beta_t-beta_expand;
% index(part1)=(beta_b-beta_t)./tmp1(part1);
index(part1)=1;
part21=logical((beta_expand<beta_t_chord).*part2);
part23=logical((beta_expand>beta_b_chord).*part2);

tmp3=beta_t-beta_chord;
tmp4=beta_chord-beta_b;
N21=(tmp1+1)./delt_beta;
N21_chord=(tmp3+1)./delt_beta;
N23=(tmp2+1)./delt_beta;
N23_chord=(tmp4+1)./delt_beta;

part22=logical(part2.*(part2-part21-part23));
N1=(beta_t_chord-beta_b+1)./delt_beta;
N3=(beta_t-beta_b_chord+1)./delt_beta;
N=N1+N3;
% index21=N21./N1;
% index23=N23_chord./N1;
index21=N21_chord./N3;
index23=N23./N3;
% index21=(N21_chord./N3+N21./N1)/2;
% index23=(N23./N3+N23./N3)/2;
index(part22)=1;
index(part21)=index21(part21);
index(part23)=index23(part23);
index=index*delt_beta;