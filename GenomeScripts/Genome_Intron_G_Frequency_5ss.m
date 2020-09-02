%SCRIPT
%To run it change the directory of fasta file; 
%Next: EDITOR/Run
blocksize=2000;
count=0;
thresh=60000;
threshNucleotide=100;
minLen=25;
begin=1;
numBlock=103;
maxSeq=min([blocksize*numBlock, 207000]);
freqC=nan(maxSeq,1);
len=nan(maxSeq,1);
nucleotideG=zeros(1, threshNucleotide);
nucleotideTotal=zeros(1, threshNucleotide);
finished = false;
tic;
while ~finished && count<maxSeq
    %read a block of sequence
    disp(['block: ' num2str(count/blocksize+1) '. There are ' num2str(numBlock-count/blocksize-1) ' blocks to go']);
    [h,s]=fastaread('C:\Users\Desktop\Genome Scripts\all_introns.fa', 'Blockread', [count+1, count+blocksize]);
    %calculate the G frequency by nucleotide
    sMat=char(nan(numel(s), threshNucleotide));
    for i=1:numel(s)
        ss=upper(s{i});
        if numel(ss)<minLen
            continue;
        end
        len(count+i)=numel(ss);        
        if numel(ss)>thresh
            freqG(count+i)=numel(find(ss(begin:thresh)=='G'))/(thresh-begin+1);
        else
            freqG(count+i)=numel(find(ss(begin:end)=='G'))/(numel(ss)-begin+1);
        end
        if numel(ss)>threshNucleotide
            sMat(i,:)=ss(1:threshNucleotide);
        else
            sMat(i,1:numel(ss))=ss;
        end
    end

    for i=1:threshNucleotide
        nucleotideG(i)=nucleotideG(i)+numel(find(sMat(:,i)=='G'));
        nucleotideTotal(i)=nucleotideTotal(i)+numel(find(sMat(:,i)));
    end
    %calculate the G freq in the first "thresh" nucleotide starting from
    %beginning
    count=count+blocksize;
    if numel(h)<blocksize
        finished=true;
    end
end
toc;
freqG=freqG(~isnan(freqG));
len=len(~isnan(len));
nucleotideG=nucleotideG./nucleotideTotal;
figure;
plot(1:threshNucleotide, nucleotideG);
title('G_{freq} by nucleotide');
xlabel('Nt position from 5SS');
ylim([0 1]);

