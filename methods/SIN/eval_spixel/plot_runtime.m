input_path = './runtime/';
suffix = '_mean_time.txt';

% semilogy seal runtime
seal_file = [input_path 'seal' suffix];
seal_data = importdata(seal_file);

ssn_file = [input_path 'ssn' suffix];
ssn_data = importdata(ssn_file);

fcn_file = [input_path 'fcn' suffix];
fcn_data = importdata(fcn_file);

our_file = [input_path 'deconv' suffix]
our_data = importdata(our_file);



% semilogy runtime
h = zeros(1,1); 
figure;
h(1)=semilogy(seal_data(:, 1),seal_data(:, 2), 'm-d','MarkerFaceColor','m', 'MarkerSize',8); hold on; 
h(2)=semilogy(ssn_data(:,1),ssn_data(:,2), 'c-s', 'MarkerFaceColor','c','MarkerSize',8); hold on;
h(3)=semilogy(fcn_data(:,1), fcn_data(:,2), 'k -^', 'MarkerFaceColor','k','MarkerSize',8); hold on;
h(4)=semilogy(our_data(:,1),our_data(:,2), 'r-o','MarkerFaceColor','r', 'MarkerSize',8); hold on;
% h(4)=semilogy(LSC_save(:,1),LSC_save(:,2), 'r-o', 'MarkerFaceColor','r','MarkerSize',8); hold on; 
% h(5)=semilogy(ERS_save(:,1),ERS_save(:,2), 'r-+', 'MarkerSize',8); hold on; 
% h(2)=semilogy(etps_save(:,1),etps_save(:,2), 'g-x', 'MarkerSize',8); hold on; 
% h(7)=semilogy(SEAL_save(:,1),SEAL_save(:,2), 'm-d','MarkerFaceColor','m', 'MarkerSize',8); hold on; 
% h(8)=semilogy(SSN_save(:,1),SSN_save(:,2), 'c-s','MarkerFaceColor','c', 'MarkerSize',8); hold on; 

hold off;
% lg = legend(h,'Ours','SLIC','SNIC','LSC','ERS', 'ETPS',  'SEAL',  'SSN', ... 
%           'Location','southeast');
lg = legend(h,'SEAL', 'SSN', 'FCN', 'Ours','Location','southeast');
set(lg,'FontSize',12, 'FontName', 'ArialMT');
xlim([200,1200])
ylim([0.005,1.5])
grid on
yticks([ 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5])
xlabel('Number of Superpixels')
ylabel('Avg. Time(log sec.)')
set(gca,'FontSize',14,'FontName', 'ArialMT')


