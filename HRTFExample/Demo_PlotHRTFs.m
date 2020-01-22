% Code to plot HRTFs
hrtf_folder = '/Volumes/MacMiniBackup/OldSSDContents/AuditoryAnalysis/golay_c57/';
%hrtf_file = '20180820_1000samp_ILDITDmod_correct_noitd.mat';
hrtf_file = '2018-08-20_1000samp_ILDITDmod_correct_noitd_ex_ILD.mat';
%hrtf_file = '2018-08-15_1000samp_ILDITDmod_correct_noitd_ex_ILD.mat';
%hrtf_file = '2018-08-16_1000samp_ILDITDmod_correct_noitd_ex_ILD.mat';

load([hrtf_folder hrtf_file], 'levels', 'peaktimes', 'allimps_percep_ILDcorrected');
api = allimps_percep_ILDcorrected;
mindb = inf;
maxdb = -inf;

titles = {'Contralateral HRTF gain', 'Ipsilateral HRTF gain', 'Contra - Ipsi'};

ahw = 5;
ahl = 7;
alw = 1;

numtone = 48;
candidate_tones = logspace(log10(5000), log10(80000), numtone * 2 + 1);
tone_edges = candidate_tones(1:2:end);
tone_centers = candidate_tones(2:2:end);


ibands = [13, 28; 29, 39; 40, 44; 45, 48];
corres = [(((0:48) - 12)')/36 tone_edges' / 1000];
corrc = {'m', 'c', 'm', 'c'};
horlines = ibands(1:3, 2);
clf


dbffts = {};
for i = 1:3
    ax(i) = subplot(1, 3, i);
    if i < 3
        irall{i} = squeeze(api(3-i, 1, :, :))';
        % this api variable has all HRIRs.
        % The structure is
        % api(ear, elevation, azimuth, time)
        % ear is 1: left, 2: right
        % elevation is 0:10:90
        % azimuth is 0:1.8:360;
        % currently, it picks up all azimuth in elevation 1 (0 degree).
        % therefore, horizontal plane HRTF.
        
        % This line does Fourier transform to get HRTFs from HRIRs.
        [freq, fftall{i}] = fftstrength(irall{i}, 500000);
        dbfft = 20 * log10(fftall{i}');
        dbfft = dbfft + 66; % to cancel arbitrary correction
        mindb = min(min(min(dbfft(1:100, 21:320))), mindb);
        maxdb = max(max(max(dbfft(1:100, 21:320))), maxdb);
        %dbfft = dbfft + 121.7719;
        dbffts{i} = dbfft;
        pcolor((0:199)* 1.8, freq / 1000, dbfft');
    else i == 3
        pcolor((0:199)* 1.8, freq / 1000, (dbffts{1} - dbffts{2})');
    end
    set(gca, 'YScale', 'log');
    colorbar
    shading flat
    colormap bluered
    ylim([10000, 80000] / 1000);
    %yticks([10, 20, 40, 60, 80])
    xlim([0, 170]);
    xlabel('RF azimuth (°)');
    title([titles{i} newline], 'FontWeight', 'normal');
    set(gca, 'TickDir', 'out');
    box off


    hold on
    for j = 1:length(horlines)
        plot(xlim, tone_centers(horlines(j)) / 1000 * [1, 1], '--', 'linewidth', 1, 'Color', 0.3 *[1 1 1]);
    end


    yticks([5, 10, 20, 30, 40, 60, 80]);
    xticks([0 50 100 150]);
    ylabel('Frequency (kHz)');

    textaxis('(dB)', [1.11, 1.1], 'fontsize', 8);
    %textaxis(texts(i), [-0.37, 1.12], 'fontsize', 12, 'fontweight', 'bold');
    for j = 1:4
        arrowaxis([-0.02, corres(ibands(j, 1))], [-0.02, corres(ibands(j, 2))], 'Color', corrc{j}, 'HeadWidth', 0, 'HeadLength', 0, 'linewidth', 3);
    end
    set(gca, 'Layer', 'top');
    set(gca, 'TickDir', 'in');

end


for i = 1:2
    set(ax(i), 'clim', [-50 50]);
end
set(ax(3), 'clim', [-60, 60]);
