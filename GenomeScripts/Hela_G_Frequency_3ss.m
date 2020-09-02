%SCRIPT
%To run it change the directory of fasta file and the csv file; 
%Next: EDITOR/Run
[h,s]=fastaread('C:\Users\Desktop\Genome Scripts\talhouarne_helaIntron.fa');
indLariat=csvread('C:\Users\Desktop\Genome Scripts\lariat_cell_hela.csv',1,0');
indh=zeros(numel(h),1);
for i=1:numel(h)
    indh(i)=sscanf(h{i}, 'hela_%d');
end
ind=zeros(numel(indLariat),1);
for i=1:numel(indLariat)
    ind(i)=find(indh'==indLariat(i));
end
sLariat=(s(ind));
sLariat=reverse(s(ind));
nLariat=numel(sLariat);
thresh=600;
begin=1;
freq=zeros(nLariat,1);
for i=1:nLariat
    ss=upper(sLariat{i});
    
    if numel(ss)>thresh
        freq(i)=numel(find(ss(begin:thresh)=='G'))/(thresh-begin+1);
    else
        freq(i)=numel(find(ss(begin:end)=='G'))/(numel(ss)-begin+1);
    end
end


%%nucleotide resolution
sLariat=reverse(s(ind));
thresh=100;
nucleotideFreq=zeros(thresh,1);
ss=nan(nLariat, thresh);
for i=1:nLariat
    if numel(sLariat{i})>thresh
        ss(i,:)=sLariat{i}(1:thresh);
    else
        ss(i,1:numel(sLariat{i}))=sLariat{i};
    end
end
for i=1:thresh
    nucleotideFreq(i)=numel(find(ss(:,i)=='G'))/numel(find(~isnan(ss(:,i))));
end

% % freq=freq(~isnan(freq));
% % len=len(~isnan(len));
% % nucleotideG=nucleotideG./nucleotideTotal;
figure;
% subplot(1,2,1);
% histogram(freq, 0:0.05:1)
% title(['G_{Freq} full length' num2str(thresh) ' nt']);
% xlabel('G_{Freq}');
% subplot(1,2,2);
plot(1:thresh, nucleotideFreq);
title('G freaquancy by nucleotide');
xlabel('Nt position from 3ss');
ylim([0 1]);
