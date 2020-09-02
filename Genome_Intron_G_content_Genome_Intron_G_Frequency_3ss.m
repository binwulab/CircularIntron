%SCRIPT
%To run it change the directory of fasta file;
%Next: EDITOR/Run
blocksize=2000;% number of sequences in one block that will be run at ones
count=0;
thresh=60000;% sequence size threshold 
threshNucleotide=100;% number of nt from 3'ss that will be ploted for G frequency
minLen=55; % minimal length of intron which will be trimmed by 30nt
begin=30;%trimming of the sequence from 3'ss side
numBlock=103;% number of total blocks that will be read
maxSeq=min([blocksize*numBlock, 207000]);
freqA=nan(maxSeq,1);
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
        ss=reverse(ss);
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
%     nucleotideG=nucleotideG+sum(sMat=='G',1);
%     nucleotideTotal=nucleotideTotal+sum(sMat~=0,1);
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
subplot(1,2,1);
hist(freqG,0.025:0.1:1);
title(['Histogram of G content']);
xlabel('G content');
subplot(1,2,2);
plot(1:threshNucleotide, nucleotideG);
title('G_{freq} by nucleotide');
xlabel('Nt position from 3SS');
ylim([0 1]);
