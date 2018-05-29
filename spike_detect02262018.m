function spike_detect02262018(input, output)

  input = strcat(input,'/**/*.h5');
  holder = dir(input);
  if exist(output,'dir')==7
  else
    mkdir(output);
  end
 addpath(genpath('/home/luly/toolbox/wave_clus-testing'),'-end');
 
  num_files = size(holder,1)
  gen_batch = cell(num_files,3);
  
  for f = 1:num_files
    base_string = holder(f).name;
    base_string = base_string(1:end-3);
    gen_batch{f,1} = strcat(holder(f).folder,'/',holder(f).name);
    gen_batch{f,2} = strcat(output,'/',base_string,'_5.mat');
    gen_batch{f,3} = holder(f).bytes;
  end
  
  for i = 1:size(gen_batch,1)
      
    if exist(gen_batch{i,2},'file')==2
      continue;
    end

%    if gen_batch{i,3} > 800000000
%      disp('Wrong file:')
%      gen_batch{i,1}
%      continue;
%    end

    file_name = gen_batch{i,1}

    all_data = h5read(file_name,'/Data/Recording_0/AnalogStream/Stream_0/ChannelData'); % raw voltage 
    info = h5read(file_name,'/Data/Recording_0/AnalogStream/Stream_0/InfoChannel'); % channel ID
    
    j = 1;
    spike_info = struct;
    conv = double(info.ConversionFactor(1)) * 10^(double(info.Exponent(1)) + 6);
    data_len = size(all_data,1)
    sr = double(1e6/info.Tick(1))
    segments_sec = double(data_len/sr);
    segments_length = double(segments_sec/60)
    param.sr = sr;
    param.segments_length = segments_length;
    param.stdmin = 5;
    param.detect_order = 4;
   
    for channel = 1:size(all_data,2)
      data = all_data(:,channel)';
      save('channel_data.mat','data');
      Get_spikes('channel_data.mat','par',param);
      load('channel_data_spikes.mat');
      if size(index,2) >= 30 
	spike_info(j).chID = info.Label(channel);
	spike_info(j).time = index';
        spike_info(j).ISI = diff(spike_info(j).time);
	spike_info(j).spikes = spikes * conv;
        spike_info(j).dV = diff(spike_info(j).spikes,1,2);
	spike_info(j).threshold = threshold * conv;
	spike_info(j).par = par;
	spike_info(j).duration = segments_sec;
	j = j+1;
      end
    end
    
    save(gen_batch{i,2},'spike_info');

  end
  
end
