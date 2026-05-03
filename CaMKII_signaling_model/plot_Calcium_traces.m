load('CA3_BK_results.mat');
% BK_con_vals=BK_con_vals([1,2,8]);
% 'r', 'b', 'g'
% colors=[1 0 0; 0 0 1; 0 1 0];
colors=jet(length(BK_con_vals));
figure;
for i=1:length(BK_con_vals)
    plot(t_gatherer{i}, (10^6)*Ca_gatherer{i}, 'LineWidth', 2, 'Color', colors(i,:));
    hold on;
end
% plot a horizontal line at 50 nm (dotted black line)
yline(50, 'k--', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Calcium Concentration (nM)');
legend([arrayfun(@(x) sprintf('BK con = %.4f mho/cm^2', x), BK_con_vals, 'UniformOutput', false)], 'Location', 'Best');
Ca_gatherer_transform = cell(size(Ca_gatherer));
% Standard deviation estimation for calcium concentration traces
Std_gatherer = zeros(size(Ca_gatherer));
FFT_gatherer= cell(1,length(BK_con_vals));
figure;
for i=1:length(BK_con_vals)
    Ca_temp=(10^6)*Ca_gatherer{i};
    % find idx b/w 2000 ms and 8000 ms
    idx=find(t_gatherer{i}>=2000 & t_gatherer{i}<=8000);
    Ca_temp=Ca_temp(idx);
    % calculate std of Ca_temp
    Ca_std=std(Ca_temp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t_unique,ia]=unique(t_gatherer{i}(idx));
    Ca_unique=Ca_temp(ia);
    t_axis_regular=2000:1:8000;
    Ca_uniform=interp1(t_unique,Ca_unique,t_axis_regular,'linear');
    fft_temp=fft(Ca_uniform);
    Fs=1000/(t_axis_regular(2)-t_axis_regular(1));
    N=length(Ca_uniform);
    f=(0:N-1)*(Fs/N);
    FFT_gatherer{i}=abs(fft_temp);
    plot(f,log10(FFT_gatherer{i}), 'LineWidth', 2, 'Color', colors(i,:));
    hold on;
    Std_gatherer(i)=Ca_std;
end
xlabel('Frequency (Hz)');
ylabel('Magnitude of FFT (log scale)');
legend([arrayfun(@(x) sprintf('BK con = %.4f mho/cm^2', x), BK_con_vals, 'UniformOutput', false)], 'Location', 'Best');
figure;
plot(BK_con_vals, Std_gatherer, 'o-', 'LineWidth', 2);
xlabel('BK con (mho/cm^2)');
ylabel('Standard Deviation of Calcium Concentration (nM)');
%%%%%% Transformation %%%%%%%%%%%
for i=1:length(BK_con_vals)
     Ca_gatherer_transform{i}=(((10^6)*Ca_gatherer{i}).^(2.72))/(5800);
end
% Plot transformed calcium concentration
figure;
for i=1:length(BK_con_vals)
    plot(t_gatherer{i}, Ca_gatherer_transform{i}, 'LineWidth', 2, 'Color', colors(i,:));
    hold on;
end
xlabel('Time (ms)');    
ylabel('Transformed Calcium Concentration');
legend([arrayfun(@(x) sprintf('BK con = %.4f mho/cm^2', x), BK_con_vals, 'UniformOutput', false)], 'Location', 'Best');
% Save the transformed calcium concentration
save('transformed_Ca_gatherer.mat', 'Ca_gatherer_transform', 't_gatherer', 'BK_con_vals');
% Plotting voltage traces: v_gatherer
figure;
for i=1:length(BK_con_vals)
    plot(t_gatherer{i}, v_gatherer{i}, 'LineWidth', 1, 'Color', colors(i,:));
    hold on;
end
xlabel('Time (ms)');    
ylabel('Voltage (mV)');
legend([arrayfun(@(x) sprintf('BK con = %.4f mho/cm^2', x), BK_con_vals, 'UniformOutput', false)], 'Location', 'Best');
