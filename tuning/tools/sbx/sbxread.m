function x = sbxread(fname,k,N,varargin)

% img = sbxread(fname,k,N,varargin)
%
% Reads from frame k to k+N-1 in file fname
% 
% fname - the file name (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.
%
% If N>1 it returns a 4D array of size = [#pmt rows cols N] 
% If N=1 it returns a 3D array of size = [#pmt rows cols]
%
% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
%
% The function also creates a global 'info' variable with additional
% informationi about the file

global info_loaded info

% check if already loaded...

if(isempty(info_loaded) || ~strcmp(fname,info_loaded))
    
    if(~isempty(info_loaded))   % try closing previous...
        try
            fclose(info.fid);
        catch
        end
    end
    
    % to be backward compatible with previous file structure...
    
%     [L,~] = strsplit(pwd,filesep);
%     if strcmp(L{end},fname)==0
%         fname = [fname filesep fname];
%     end

    try
        load(fname);
    catch
        error('No such file %s', fname)
    end
    
    if(exist([fname ,'.xevents'],'file')) % is there some external ttl event data
        xefid = fopen([fname ,'.xevents'],'r');
        edata = uint8(fread(xefid,inf,'uint8'));
        info.xevents.code = edata;
        info.xevents.binary = dec2bin(edata,8);
    else
        info.xevents.code = [];
        info.xevents.binary = [];
    end

    if(exist([fname ,'.align'],'file')) % aligned?
        info.aligned = load([fname ,'.align'],'-mat');
    else
        info.aligned = [];
    end   
    
    info_loaded = fname;
    
    if(~isfield(info,'sz'))
        info.sz = [512 796];    % it was only sz = .... 
    end
    
    if(~isfield(info,'scanmode'))
        info.scanmode = 1;      % unidirectional
    end
    
    if(info.scanmode==0)
        info.recordsPerBuffer = info.recordsPerBuffer*2;
    end
    
    if(~isfield(info,'chan'))
        switch info.channels
            case 1
                info.nchan = 2;      % both PMT0 & 1
                factor = 1;
            case 2
                info.nchan = 1;      % PMT 0
                factor = 2;
            case 3
                info.nchan = 1;      % PMT 1
                factor = 2;
        end
    else
        info.nchan = info.chan.nchan;
    end
    
    info.fid = fopen([fname '.sbx']);
    d = dir([fname '.sbx']);
    
    info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record 
    
    if isfield(info,'scanbox_version') 
        switch(info.scanbox_version)
            case 2
            info.max_idx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 - 1;
            info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record 
            case 3
            info.max_idx = d.bytes/prod(info.sz)/info.nchan/2 -1;
            info.nsamples = prod(info.sz)*info.nchan*2;
            otherwise
                error('Invalid Scanbox version');
        end
    else
        info.max_idx =  d.bytes/info.bytesPerBuffer*factor - 1;
    end
end

if(isfield(info,'fid') && info.fid ~= -1)
    
    % nsamples = info.postTriggerSamples * info.recordsPerBuffer;
        
    try
        fseek(info.fid,k*info.nsamples,'bof');
        x = fread(info.fid,info.nsamples/2 * N,'uint16=>uint16');
        x = reshape(x,[info.nchan info.sz(2) info.recordsPerBuffer  N]);
    catch
        error('Cannot read frame.  Index range likely outside of bounds.');
    end

    x = intmax('uint16')-permute(x,[1 3 2 4]);
    
    % folding lines
    
    if isfield(info,'fold_lines')
        if info.fold_lines>0
            sx = size(x);
            if length(sx) == 3
                sx(4)=1;
            end
            x = reshape(x,sx(1),info.fold_lines,[],sx(3),sx(4));
            x = permute(x,[1 2 4 3 5]);
            sx = size(x);
            x = reshape(x,sx(1),sx(2),sx(3),[]);
            
% %             t = info.frame(idx)*512+info.line(idx);   % recode TTLs
% %             info.frame = floor(t/info.fold_lines);
% %             info.line = mod(t,info.fold_lines);
            
        end
               
    end
    
else
    x = [];
end