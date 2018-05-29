function  PS_hist_lin(varargin)

  pvpmod(varargin);

  if ~exist('input', 'var')
    input = '.';
  end

  input = strcat(input,'/*.mat');

  holder = dir(input);
  num_files = size(holder,1);
  for f = 1:num_files
    base_string = holder(f).name
    base_string = base_string(1:end-4);
    csvfile = strcat(base_string,'_PS_Cpeak.csv');

    if exist(csvfile,'file')
      continue;
    end

    pdffile = strcat(base_string,'_PS_hist_Cpeak.pdf');

    if exist(pdffile,'file')
      continue;
    end

    file_name = strcat(holder(f).folder,'/',holder(f).name)
    load(file_name);

    if size(spike_info,2) ==1
      continue;
    end
    duration = spike_info(1).duration
    Crosscorr = struct;

    k = 1;
    pdff = {};

    for i = 1:size(spike_info,2)
      if length(spike_info(i).time) < 10
        continue;
      end

      for j = 1:size(spike_info,2)
        if length(spike_info(j).time) < 10
          continue;
        end
         if i == j
          continue;
         end
        k1 = char(spike_info(i).chID{1});
        k2 = char(spike_info(j).chID{1});
        h1 = histogram(spike_info(i).time,0:1000:duration*1000);
        mean_fire_rate_ref = mean(h1.Values);
        h2 = histogram(spike_info(j).time,0:1000:duration*1000);
        mean_fire_rate_fow = mean(h2.Values);

        xalph1 = double(k1(1));
        xalph2 = double(k2(1));
        if xalph1 < 73
          x1 = double(k1(1)-64);
        else
          x1 = double(k1(1)-65);
        end
        if xalph2 < 73
          x2 = double(k2(1)-64);
        else
          x2 = double(k2(1)-65);
        end
        y1 = str2double(k1(2:length(k1)));
        y2 = str2double(k2(2:length(k2)));
        d_1 = abs((x2-x1) * 100);
        d_2 = abs((y2-y1) * 100);
        d = abs(sqrt(d_1.^2 + d_2.^2));
        n1 = length(spike_info(i).time);
        if n1 < 30
            continue;
        end
        n2 = length(spike_info(j).time);
        if n2 < 30
            continue;
        end
        
        [Q0, deltaT] = CCC(spike_info(i).time,spike_info(j).time,1,50,duration,0);
         [height,pos] = max(Q0);
        pos = deltaT(pos);
        height = round(height);
        if isempty(Q0) | height < 80
          clear Q0;
          clear deltaT;
          continue;
        end
        [deltaT3,Q3,ic] = unique(deltaT);
        Q3 = Q0(Q3);
        
        %PS
        deltaT2 = deltaT;
        deltaT2(deltaT >15 | deltaT < -15) = [];
        
        deltaT(deltaT >2 | deltaT < -2) = [];
        
        psht = histogram(deltaT, -2:0.1:2,'Normalization','probability');
        ps1 = psht.BinCounts;
        [ps_event,PS_n] = max(ps1);
        C_probability = round(psht.Values(PS_n),3);
        if C_probability < 0.3
            continue;
        end
        C_peak =round(ps_event/sqrt(n1*n2),3);
        
        PS_deltaT = mode(deltaT);
        if PS_deltaT < 0 
            continue;
        end
        ps_r = round(ps_event/n1,3);
        
        subplot(3,1,1);
        mdistanc =  max(deltaT3) - min(deltaT3) -0.05;
        if mdistanc > 6
          mdistanc = 6;
        elseif mdistanc < 0
          continue;
        end
        if length(deltaT3) < 3
            continue;
        end
        findpeaks(Q3,double(deltaT3),'MinPeakHeight',75,'MinPeakProminence',60,'MinPeakDistance',mdistanc,'Annotate','extents');
        clear Q3;
        clear Q0;
        clear deltaT3;
        title([num2str(height),'|',num2str(pos)]);
        subplot(3,1,2);
        histogram(deltaT2, -15:1:15,'Normalization','probability');
        clear deltaT2;
        title([k1,' <-> ',k2,num2str(C_probability),'|',num2str(C_peak)]);
        subplot(3,1,3);
        histogram(deltaT, -2:0.1:2,'Normalization','probability');
        clear deltaT;
        title([k1,' <-> ',k2,' | ',num2str(PS_deltaT),' | ',num2str(ps_event),' | ',num2str(C_probability) ' | Ratio=',num2str(ps_r)]);

        pdff{k} = strcat('tmp_',num2str(k),'.pdf');
        print(pdff{k}(1:end-4),'-dpdf','-fillpage');
        Crosscorr(k).id = k;
        Crosscorr(k).chID1 = k1;
        Crosscorr(k).chID2 = k2;
        Crosscorr(k).Dist = d;
        Crosscorr(k).chID1_spikes = n1;
        Crosscorr(k).chID2_spikes = n2;
      %[h,p1] = kstest((norm_offsets)');                                 
 
        Crosscorr(k).PS_event = ps_event;
        Crosscorr(k).PS_DeltaT = PS_deltaT;
        Crosscorr(k).C_probability = C_probability;
        Crosscorr(k).C_peak = C_peak;
        Crosscorr(k).PS_HZ = round((ps_event)/180,3);
        Crosscorr(k).firing_ref = round(n1/180,3);
        Crosscorr(k).firing_fow = round(n2/180,3);
        Crosscorr(k).PS_others = sum(ps1)-ps_event;
        Crosscorr(k).PS_Ratio = ps_r;
        Crosscorr(k).ref_ISI = mean(spike_info(i).ISI);
        Crosscorr(k).tag_ISI = mean(spike_info(j).ISI);
        Crosscorr(k).mean_fire_rate_ref = mean_fire_rate_ref;
        Crosscorr(k).mean_fire_rate_fow = mean_fire_rate_fow;
         k = k+1;
      end
    end
    if k > 1  % or if size(Crosscorr,2) > 0                                                        
      B = struct2table(Crosscorr);
      [B,ind] = sortrows(B, 'PS_event', 'descend');
      ID = zeros(size(B,1),1);
      ID(:,1) = 1:size(B,1);
      B.id=ID;
      writetable(B, csvfile);
      append_pdfs(pdffile,pdff{ind});
      delete tmp*.pdf;
      clear B;
    end
    clear Crosscorr;
    clear pdff;

 end



    
            
        
       
       
