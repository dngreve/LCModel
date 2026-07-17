classdef lcmodel < handle
  properties
    basis_file = '';
    outdir = '';
    met = [];
    ref = [];
    Fs = 0; % sample frequency
    ppmst  = 4.0; % "start" freq in ppm
    ppmend = 1.9; % "end" freq in ppm
    nunfil = 0; % Nt = size(fid_met,1);
    hzpppm = 0; % 123.137;
    echot = 0; % echo time 17.6
    %deltat = 1/Fs
    DoWS = 'T'; % Do Water Scaling
    DoECC = 'F'; % Do eddy current correction
    sddegp=0; % expected std of degppm
    nomit=0; % n-omit - n metabolites to exclude
    dorefs1='F'; % 'DOREFS(1)=F' T=Use water peak as landmark
    atth2o = 1; % scale used for water scaling
    wconc = 39590; % nmr-visible water concentration (mM); tissuetype dep
    id = '29'; % optional string for identifying in plots
    fmtdat = '2E15.6'; % itself must be in ''
    %fprintf(fp2, '%15.6e%15.6e\n',interleaveRealImag(fid_met));

    % These are related to multi-voxel
    ndcols = 0;
    ndrows = 0;
    ndslic = 0;
    icolst = 0;
    icolen = 0;
    irowst = 0;
    irowen = 0;
    islice = 0;
    nvoxsk= 0

    % volume and tramp are used to tramp/volume
    volume = '2.700e+01'; % VOI/voxel volume
    tramp = 1; % inv related to gain

    title = 'title';
    projectname = ''; % prepend
    lcsv = 11; %11= make .csv file, 0 not

    % Probably best to set these next four to 0 for voxel-wise
    % analysis to reduce the number of files created
    ltable = 7; % 7=make the .table file, 0 not
    lps = 8; %8=make .ps file, 0 not
    lprint = 6; %6= make .print file, 0 not
    lcoord = 9; %9=make .coord file, 0 not
    pgnorm = sprintf('%cUS%c',char(39),char(39)); % page style
    ipage2 = 1; %0=suppress printing of a 2nd page
    neach = 99; % Number of metabolites for which individual plots are to be made.

    filpri = '';
    filcoo = '';
    bash_script = '';
    srcraw = '';
    savedir = '';
    key=210387309;
    enableLCMspecFileOutput = 1;
    enableLCMdetailedParametersOutput = 1;
    enableLCMallSpecOutput = 1;
    crs0 = []; % not part of the control file
  end % properties
  
  methods

    function lcmfid = midas2lcmodel(lcm,midasspectrum,fidtype)
      % Convert midas specra to lcmodel fid. 
      % Insert midas spectrum into a two sided spectrum. For met,
      % the other side is 0s. For ref, the midas spectrum is
      % centered with 0s on either size
      [nfid nvox] = size(midasspectrum);
      spec2 = complex(zeros(2*nfid,nvox),zeros(2*nfid,nvox));
      ind = 1:nfid; % ref has a different center
      if(strcmp(fidtype,'ref')) ind = round(nfid/2)+ind; end
      spec2(ind,:) = midasspectrum;
      % Compute the FID
      lcmfid = ifft(spec2);
      lcmfid = lcmfid(1:nfid,:);
    end

    function write_fid(lcm,fidtype)
      % Note: can't have more than 65536 voxels in the file
      mkdirp(lcm.outdir);
      outfile = sprintf('%s/%s.%s',lcm.outdir,lcm.projectname,fidtype);
      fp=fopen(outfile,'w');
      fprintf(fp,' $SEQPAR\n');
      fprintf(fp,' echot=%6.2f\n', lcm.echot); % echo time
      fprintf(fp,' hzpppm=%12.4e\n', lcm.hzpppm);
      fprintf(fp,' $END\n');
      fprintf(fp,' $NMID\n');
      fprintf(fp,' id= %s, fmtdat= ''(%s)''\n',lcm.id,lcm.fmtdat);
      fprintf(fp,' volume = %s\n',lcm.volume);
      fprintf(fp,' tramp= %g\n',lcm.tramp);
      fprintf(fp,' $END\n');
      if(strcmp(fidtype,'met')) fid = lcm.met; end
      if(strcmp(fidtype,'ref')) fid = lcm.ref; end
      fprintf(fp, '%15.6e%15.6e\n',[real(fid(:)) imag(fid(:))]');
      fclose(fp);
      return
    end
    
    function write_control(lcm,nsegs)
      % write control file. Also calls write_fid
      % nsegs is the number of partitions, typically the number in
      % the mask
      [lcm.nunfil lcm.ndcols] = size(lcm.met);
      nvox = lcm.ndcols;
      if(~exist('nsegs','var')) nsegs = []; end
      if(nvox > 65536 | ~isempty(nsegs))
        met0 = lcm.met;
        ref0 = lcm.ref;
        outdir0 = lcm.outdir;
        if(isempty(nsegs)) 
          seglength = 65536;
          nsegs = ceil(nvox/seglength); 
        else
          seglength = floor(nvox/nsegs);
        end
        fprintf('nsegs = %d seglength = %d nvox = %d\n',nsegs,seglength,nvox);
        fname = sprintf('%s/vox.txt',outdir0);
        fp = fopen(fname,'w');
        fprintf(fp,'%d %d %d\n',nsegs,seglength,nvox);
        fclose(fp);
	nprint = round(nsegs/10);
        tic;
        for nthseg = 1:nsegs
          % the v%06 folder will be the index within the mask
          % starting at 1, so goes from 1 to nmask. 
          lcm.outdir = sprintf('%s/v%06d',outdir0,nthseg);
          controlfile = sprintf('%s/%s.control',lcm.outdir,lcm.projectname);
          %if(fast_fileexists(controlfile)) continue; end
          mkdirp(lcm.outdir);
          n1 = (nthseg-1)*seglength + 1;
          n2 = (nthseg)*seglength;
          if(n2 > nvox) n2 = nvox; end
          if(nthseg==1|rem(nthseg,nprint)==0) fprintf('seg %d %d %d  %6.3f\n',nthseg,n1,n2,toc/60); end
          fname = sprintf('%s/vox.indices.dat',lcm.outdir);
          fp = fopen(fname,'w');
          fprintf(fp,'%d %d\n',n1,n2);
          fclose(fp);
          lcm.met = met0(:,n1:n2);
          lcm.ref = ref0(:,n1:n2);
          lcm.write_control();
        end
        lcm.outdir = outdir0;
        lcm.met = met0;
        lcm.ref = ref0;
        return
      end
      % Does not get here if nvox>65536 
      % check ref dim = met dim
      mkdirp(lcm.outdir);
      controlfile = sprintf('%s/%s.control',lcm.outdir,lcm.projectname);
      fp=fopen(controlfile,'w');
      fprintf(fp,' $LCMODL\n');
      fprintf(fp,' key=%d\n',lcm.key);
      fprintf(fp,' title=''%s''\n', lcm.projectname );
      fprintf(fp,' srcraw=''%s''\n', lcm.srcraw );
      fprintf(fp,' savdir=''%s''\n', lcm.outdir); %?
      % This setting may cause an init error, eg, if it is too wide
      fprintf(fp,' ppmst=%f\n', lcm.ppmst);  
      fprintf(fp,' ppmend=%f\n', lcm.ppmend);

      fprintf(fp,' nunfil = %d\n', lcm.nunfil); 
      fprintf(fp,' ltable = %d\n',lcm.ltable);
      fprintf(fp,' lps = %d\n',lcm.lps);
      fprintf(fp,' lcsv = %d\n',lcm.lcsv);
      fprintf(fp,' hzpppm = %12.4e\n', lcm.hzpppm); 
      fprintf(fp,' filtab = ''%s/%s.table''\n',lcm.outdir,lcm.projectname);
      fprintf(fp,' filraw = ''%s/%s.met''\n',lcm.outdir,lcm.projectname);
      fprintf(fp,' filh2o = ''%s/%s.ref''\n',lcm.outdir,lcm.projectname);
      fprintf(fp,' filps = ''%s/%s.ps''\n',lcm.outdir,lcm.projectname);
      fprintf(fp,' filcsv = ''%s/%s.csv''\n',lcm.outdir,lcm.projectname);
      fprintf(fp,' filbas=''%s''\n', lcm.basis_file);
      if lcm.enableLCMdetailedParametersOutput
        fprintf(fp,' filpri = ''%s/%s.log''\n',lcm.outdir,lcm.projectname);
        fprintf(fp,' lprint = %d\n',lcm.lprint);
      end
      if lcm.enableLCMspecFileOutput
        fprintf(fp,' filcoo = ''%s/%s.spec''\n',lcm.outdir,lcm.projectname);
        fprintf(fp,' lcoord= %d\n',lcm.lcoord);
      end
      if lcm.enableLCMallSpecOutput
        fprintf(fp,' neach = %d\n',lcm.neach);
      end
      fprintf(fp,' echot=%6.2f\n', lcm.echot);
      fprintf(fp,' deltat=%11.3e\n', 1/lcm.Fs);
      fprintf(fp,' dows = %s\n',lcm.DoWS);
      fprintf(fp,' doecc = %s\n',lcm.DoECC);
      fprintf(fp,' SDDEGP=%g\n',lcm.sddegp);
      fprintf(fp,' NOMIT=%d\n',lcm.nomit);
      fprintf(fp,' DOREFS(1)=%s\n',lcm.dorefs1);
      fprintf(fp,' ATTH2O=%d\n',lcm.atth2o);
      fprintf(fp,' WCONC=%g\n',lcm.wconc);
      fprintf(fp,' PGNORM=%s\n',lcm.pgnorm);
      fprintf(fp,' IPAGE2=%d\n',lcm.ipage2);
      fprintf(fp,' ndcols = %d\n',lcm.ndcols);
      fprintf(fp,' icolen = %d\n',lcm.ndcols);
      fprintf(fp,' ndrows = %d\n',1);
      fprintf(fp,' ndslic = %d\n',1);
      fprintf(fp,' nvoxsk = %d\n',0);
      fprintf(fp,' $END\n');
      fclose(fp);
      lcm.write_fid('met');
      lcm.write_fid('ref');
      return;
    end

    function mri = loadresults(lcm,lcmdir)
      fname = sprintf('%s/mask.nii.gz',lcmdir);
      mask = MRIread(fname);
      indmask = find(mask.vol);
      % cd ldcmdir; cat */dng.all.dat > dng.dat
      fname = sprintf('%s/dng.all.dat',lcmdir);
      dat = load(fname);
      indmasklcm = dat(:,1);
      nmeasdat = size(dat,2);
      % Dont include the first 3: index, col, row
      bm = zeros(nmeasdat-3,prod(mask.volsize));
      bm(:,indmask(indmasklcm)) = dat(:,4:end)';
      mri = mask;
      mri.vol = fast_mat2vol(bm,mri.volsize);
    end
    
    function mri = loadresults2(lcm,lcmdir)
      field = 'NAA_Cr';
      fname = sprintf('%s/mask.nii.gz',lcmdir);
      mask = MRIread(fname);
      indmask = find(mask.vol);
      fname = sprintf('%s/vox.txt',lcmdir);
      segs = load(fname);
      nsegs = segs(1);
      mri = mask;
      mri.vol = zeros(size(mask.vol));
      tic;
      for segno = 1:nsegs
        fname = sprintf('%s/v%06d/dng.csv',lcmdir,segno);
        if(~exist(fname,'file')) continue; end
        a = readtable(fname);
        if(isempty(a)) continue; end
        fprintf('%d %g -----------------------\n',segno,toc);
        fname = sprintf('%s/v%06d/vox.indices.dat',lcmdir,segno);
        segindices = load(fname); %1-based
        VarNames = char(a.Properties.VariableNames);
        indnaacr = strmatch(field,VarNames);
        voxno = a{:,2};  %1-based
        ind = indmask(voxno+segindices(1)-1);
        mri.vol(ind) = a{:,indnaacr};
      end
    end
    
  end % methods

end





