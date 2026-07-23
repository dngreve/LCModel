function sid2lcmodel(topostr,surficostr,subjectstr)

%topo = 4; % 1=vol, 2=roi, 3=roitfc 4=lhsurf, 5=rhsurf
%surfico = 5;

% sid2lcmodel.m -- see also metasurfer2lcmodel.m
% has nothing to do with sid, should change the name
% takes result of metasurfer-pp, extracts each voxel in the mask,
% and creates a folder for lcmodel for each voxel
% next: cat */dng.dat > dng.dat


% addpath /autofs/space/iddhi_002/users/greve/projects/mrsi/lcmodel

%msdir = '/autofs/space/iddhi_005/users/greve/epsi/zip2/iwbs-001/mrsi/metasurfer';
%msdir = '/autofs/space/iddhi_005/users/greve/ku-mrsi/source/146588_EGG_018_01/metasurfer.parlite';

topdir = '/autofs/space/iddhi_005/users/greve/ku-mrsi';
msdirname = 'metasurfer.parlite.v2';

topo = sscanf(topostr,'%d');
surfico = sscanf(surficostr,'%d'); % irrelevant if not a surface
randvox = 0;

% Set mergevox=0 to have different files for each voxel. 
% Set mergevox=1 to have one file for all voxels.
% if mergevox>1, then it only keeps mergevox voxels.
mergevox = 0;

if(~exist('subjectstr','var')) subjectstr = ''; end

if(isempty(subjectstr) | strcmp(subjectstr,'all'))
  fname = sprintf('%s/subjects.txt',topdir);
  slist = char(textread(fname,'%s'));
  fprintf('Processing all subjects\n');
else
  slist = subjectstr;
  fprintf('Processing subject %s\n',slist);
end
nslist = size(slist,1);
fprintf('%s ns=%d, topo=%d surfico=%d\n',msdirname,nslist,topo,surfico);

tic
for ns = [1:nslist]
  subject = deblank(slist(ns,:));
  if(strcmp(subject,'egg-011-01')) 
    fprintf('Skipping %s\n',subject);
    continue
  end
  msdir = sprintf('%s/%s/%s',topdir,msdirname,subject);
  fprintf('%2d %s %6.4f ========================\n',ns,subject,toc/60);

  if(topo == 1) %Vol
    maskname = sprintf('%s/brainmask.nii.gz',msdir);
    metrealname =  sprintf('%s/met.real.nii.gz',msdir);
    metimagname =  sprintf('%s/met.imag.nii.gz',msdir);
    refrealname =  sprintf('%s/ref.real.nii.gz',msdir);
    % Below is a bug fixed on 7/17/26. It caused the imag ref to be equal to the
    % real ref. It affects the absolute measures (eg, NAA) but, shockingly,
    % has almost no impact on ratios (eg, NAA/Cr).  Claude explained it to
    % me, see below.
    %refimagname =  sprintf('%s/ref.real.nii.gz',msdir);
    refimagname =  sprintf('%s/ref.imag.nii.gz',msdir);
    outtop = sprintf('%s/vol.lcm',msdir);
  end
  if(topo == 2 || topo == 3) %ROI/GTM
    if(topo == 2) gtmdir = sprintf('%s/gtm/notfc',msdir); end
    if(topo == 3) gtmdir = sprintf('%s/gtm/tfc',msdir); end
    maskname = '';
    metrealname =  sprintf('%s/met.real/gtm.nii.gz',gtmdir);
    metimagname =  sprintf('%s/met.imag/gtm.nii.gz',gtmdir);
    refrealname =  sprintf('%s/ref.real/gtm.nii.gz',gtmdir);
    %refimagname =  sprintf('%s/ref.real/gtm.nii.gz',gtmdir); % bug, fixed 7/17/26
    refimagname =  sprintf('%s/ref.imag/gtm.nii.gz',gtmdir);
    outtop = sprintf('%s/lcm',gtmdir);
  end
  if(topo == 4 | topo == 5) % Surf
    if(topo == 4) hemi = 'lh'; end
    if(topo == 5) hemi = 'rh'; end
    maskname = sprintf('%s/surf%d/mask.%s.nii.gz',msdir,surfico,hemi);
    metrealname =  sprintf('%s/surf%d/met.real.%s.sm10.nii.gz',msdir,surfico,hemi);
    metimagname =  sprintf('%s/surf%d/met.imag.%s.sm10.nii.gz',msdir,surfico,hemi);
    refrealname =  sprintf('%s/surf%d/ref.real.%s.sm10.nii.gz',msdir,surfico,hemi);
    refimagname =  sprintf('%s/surf%d/ref.imag.%s.sm10.nii.gz',msdir,surfico,hemi);
    outtop = sprintf('%s/surf%d/lcm.%s.sm10',msdir,surfico,hemi);
  end
  
  if(0)
    % To read in the data and write out a volume
    lcm = lcmodel;
    a = lcm.loadresults(outtop);
    fname = sprintf('%s/lcm.nii.gz',outtop);
    MRIwrite(a,fname);
    continue
  end
  
  metreal = MRIread(metrealname);
  metimag = MRIread(metimagname);
  refreal = MRIread(refrealname);
  refimag = MRIread(refimagname);
  if(mergevox > 1)
    metreal.vol = metreal.vol(:,1:mergevox,:,:);
    metimag.vol = metimag.vol(:,1:mergevox,:,:);
    refreal.vol = refreal.vol(:,1:mergevox,:,:);
    refimag.vol = refimag.vol(:,1:mergevox,:,:);
    metreal.volsize(2) = mergevox; metreal.width = mergevox; metreal.nvoxels = mergevox;
    metimag.volsize(2) = mergevox; metimag.width = mergevox; metimag.nvoxels = mergevox;
    refreal.volsize(2) = mergevox; refreal.width = mergevox; refreal.nvoxels = mergevox;
    refimag.volsize(2) = mergevox; refimag.width = mergevox; refimag.nvoxels = mergevox;
  end
  met = metreal.vol + i*metimag.vol;
  ref = refreal.vol + i*refimag.vol;
  
  if(~isempty(maskname)) 
    mask = MRIread(maskname); 
  else
    mask = metreal;
    mask.vol = ones(metreal.volsize);
  end
  indmask = find(mask.vol);
  nmask = length(indmask);
  %[rb cb sb] = ind2sub(mask.volsize,indmask);%1-based
  [rmask1 cmask1 smask1] = ind2sub(metreal.volsize,indmask);
  crs0 = [cmask1 rmask1 smask1]-1;
  
  nspect = metreal.nframes;
  
  if(~mergevox) outtop = sprintf('%s.pv',outtop); end
  if(mergevox) outtop = sprintf('%s.mv',outtop); end
  if(randvox) 
    outtop = sprintf('%s.rand',outtop); 
    % Set it up so that the N voxels are invalid. If N>=15, then
    % lcmodel stops and does not process the remaining voxels
    nn = reshape1d([51:100; 50:-1:1]);
    met(:,nn(50+[1:2]),:,200) = 10e10;
    %met(:,nn(5+[1:15]),:,200) = 10e10;
    %met(:,[1 end],:,200) = 10e10;
  end
  mkdirp(outtop);
  outtop
  size(met)
  
  echot = 15.6; % From Echo_Time in subject.xml for AD/TDAD/EGG
                %echot = 17.6; % From Echo_Time in subject.xml for other study
  
  % Select based on echot
  basis_file = 'sim_se_csi_te16.basis';
  
  Precession_Frequency = 123.231196; % from subject.xml for AD/TDAD/EGG
                                     %Precession_Frequency =  123.136807; % other study
  
  % Real_Dwell_Time. This value is at the unit of 0.1 µs, thus 4000
  % means 400 µs. Since two echos (positive and negative echo) are
  % averaged during processing, the effective spectral width should be
  % 1/800 = 1.25 kHz
  RDT = 6300;% Value straight from the subject.xml 400
  RDT = 4000;% Value straight from the subject.xml 
  DELTAT=2*RDT*1e-7; 
  %DELTAT=0.00104058; % other study DELTAT=0.00104058 = 961Hz???
  Fs = 1/DELTAT; 
  
  %fname = 'seg0010/dng.csv';
  %fp = fopen(fname,'r');
  %s1 = fgets(fp);
  %s2 = fgets(fp);
  %fclose(fp);
  
  lcm = lcmodel;
  lcm.ltable = 0;
  lcm.lps = 0;
  lcm.lprint = 0;
  lcm.lcoord = 0;

  lcm.basis_file = sprintf('/homes/4/greve/.lcmodel/basis-sets/%s',basis_file);
  lcm.Fs = Fs;
  lcm.hzpppm = Precession_Frequency;
  lcm.echot = echot;
  lcm.projectname = 'dng';
  lcm.id = '28';
  lcm.outdir = outtop;
  lcmmet = lcm.midas2lcmodel(fast_vol2mat(met),'met'); %permute(squeeze(met),[2 1])
  lcmref = lcm.midas2lcmodel(fast_vol2mat(ref),'ref'); %permute(squeeze(ref),[2 1])
  lcm.crs0 = crs0;
  
  vlcm = metreal; % template
  vlcm.vol = real(fast_mat2vol(lcmmet,vlcm.volsize));
  fname = sprintf('%s/lcm.met.real.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  vlcm.vol = imag(fast_mat2vol(lcmmet,vlcm.volsize));
  fname = sprintf('%s/lcm.met.imag.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  vlcm.vol = abs(fast_mat2vol(lcmmet,vlcm.volsize));
  fname = sprintf('%s/lcm.met.mag.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  
  vlcm.vol = real(fast_mat2vol(lcmref,vlcm.volsize));
  fname = sprintf('%s/lcm.ref.real.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  vlcm.vol = imag(fast_mat2vol(lcmref,vlcm.volsize));
  fname = sprintf('%s/lcm.ref.imag.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  vlcm.vol = abs(fast_mat2vol(lcmref,vlcm.volsize));
  fname = sprintf('%s/lcm.ref.mag.nii.gz',lcm.outdir);
  MRIwrite(vlcm,fname);
  
  fname = sprintf('%s/mask.crs.dat',lcm.outdir);
  fp = fopen(fname,'w');
  fprintf(fp,'%3d %3d %3d\n',[cmask1' rmask1' smask1']'-1);
  fclose(fp);
  fname = sprintf('%s/nmask.dat',lcm.outdir);
  fp = fopen(fname,'w');
  fprintf(fp,'%6d\n',length(indmask));
  fclose(fp);

  fname = sprintf('%s/mask.nii.gz',lcm.outdir);
  MRIwrite(mask,fname);
  
  lcm.met = lcmmet(:,indmask);
  lcm.ref = lcmref(:,indmask);
  if(mergevox == 0) lcm.write_control(nmask);
  else              lcm.write_control(1);
  end
  
end

fprintf('sid2lcmodel done\n');

return

% This is claude's answer as to why the ratios did not change when the
% imag part of the ref was really set to the real.

% Since every metabolite's concentration is fit as a linear-combination
% coefficient against this uniformly-rescaled datat, every metabolite
% inherits the same fcalib factor — right or wrong. A ratio like NAA/Cr is
% (C_NAA·fcalib)/(C_Cr·fcalib) — fcalib cancels algebraically, exactly,
% regardless of what value it took. This is precisely why the MRS field
% leans so heavily on ratio reporting — it's structurally immune to exactly
% this class of error (water-referencing uncertainty, coil sensitivity,
% partial volume, and, as you demonstrated, a real software bug).

% Why it wasn't total garbage: area_water comes from areawa() (lines
% 5413–5487), which fits a log-linear regression to alog(cabs(h2ot(j)))
% across the water FID's time samples — it only ever looks at the magnitude
% of the complex water signal via cabs(). Your bug produced ref = refreal +
% i·refreal = (1+i)·refreal — a deterministic, smooth distortion (a fixed
% complex rotation-and-scaling, |1+i| = √2, applied identically to every
% real, physically-smooth, decaying sample), not noise or
% instability. cabs((1+i)·x) = √2·|x| is still smooth, still positive — so
% the regression still converges to a well-defined (if numerically wrong)
% area_water, which flows into a well-defined (if wrong) fcalib, which
% produces well-defined (if wrong) absolute concentrations. Nothing hit a
% .le. 0. guard or produced NaN — the bug corrupted the value, not the
% validity, of the intermediate quantities. Had the bug instead introduced
% something incoherent (random noise, a mid-FID sign flip, a discontinuity),
% you likely would have seen garbage or a crash — it's the specific shape of
% this bug (a fixed linear transform of a smooth real signal) that made it
% "quietly wrong" rather than catastrophic.


tic;
fprintf('nbrain = %d\n',nbrain);
for vno = 1:nbrain
  fprintf('vno %d %d  (%d,%d,%d) %4.1f\n',vno,indbrain(vno),cb(vno), rb(vno), sb(vno),toc/60);
  crs1 = [cb(vno) rb(vno) sb(vno)];
  %rcs1 = crs1([2 1 3]); % don't add 1
  rcs1 = [1 cb(vno) 1]; % for surfaces
  
  % In midas, the spectrum is stored with the highest freq first, so
  % flip it to make the lowest frequency first, or not 
  refspectrum0 = (squeeze(ref(rcs1(1),rcs1(2),rcs1(3),:)));
  metspectrum0 = (squeeze(met(rcs1(1),rcs1(2),rcs1(3),:)));
  
  spec2 = complex(zeros(1024,1),zeros(1024,1)); %extended spectrum to be double width
  spec2(1:512) = metspectrum0; %add 512 zeros to the end of spectrum
  fid_met = ifft(spec2);
  %fid_met(513:1024) = 0; %replace 512 points with zeros at the tail part of the FID
  fid_met = fid_met(1:512);

  %water reference
  specRef2 = complex(zeros(1024,1),zeros(1024,1));
  specRef2(257:768) = refspectrum0; %extend water spectrum both sides by 256 points
  fid_ref = ifft(specRef2);
  %fidReft(513:1024) = 0; %replace 512 points with zeros at the tail part of the FID
  fid_ref = fid_ref(1:512);
  
  outdir = sprintf('%s/%s/v%06d.%03d.%03d/',mrsidir,outtop/crs1(1),crs1(2),crs1(3));
  mkdirp(outdir);

  fname = sprintf('%s/crs1.dat',outdir);
  fp = fopen(fname,'w');
  fprintf(fp,'%3d %3d %3d  %5d\n',crs1(1),crs1(2),crs1(3),indbrain(vno));
  fclose(fp);
  
  fname = sprintf('%s/spectra.dat',outdir);
  fp = fopen(fname,'w');
  fprintf(fp,'%g %g\n',[real(metspectrum0) imag(metspectrum0)]');
  fclose(fp);
  
  projectName = 'dng';
  localAccessDir = outdir;
  selfAccessDir = outdir;
  lcmodel_file_writer(fid_met,fid_ref,Fs,projectName,basis_file,localAccessDir,selfAccessDir);

 
  fname = sprintf('%s/run.csh',outdir);
  fp = fopen(fname,'w');
  fprintf(fp,'#!/bin/csh -f\n');
  fprintf(fp,'set ctrl = %s/dng.control\n',outdir);
  fprintf(fp,'set ps = %s/dng.ps\n',outdir);
  fprintf(fp,'set ud = `UpdateNeeded $ps $ctrl`\n');
  fprintf(fp,'if(! $ud) exit 0\n');
  fprintf(fp,'echo starting lcmodel \n');
  fprintf(fp,'rm -f %s/lcm.done\n',outdir);
  fprintf(fp,'date > %s/lcm.started\n',outdir);
  fprintf(fp,'fs_time lcmodel < %s/dng.control\n',outdir);
  fprintf(fp,'date > %s/lcm.done\n',outdir);
  fprintf(fp,'echo done\n');
  fclose(fp);
end

return
