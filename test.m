System = 'data/S0_50N0_10SIM';

% Load 
load([System '/fov1_track.mat']);
% Calculate MSD
knownMSD = calcMSD(Tracks,16.25,0);
load([System '/movie4QMData.mat'])
QMMSD = CorrectedMSD;
preTracks = Tracks;
preTracks(:,3:4) = [];
preMSD = calcMSD(preTracks,16.25,0);
%Plot MSDs

fig1 = figure();
whitebg(fig1,[1,1,1])
scatter(1:size(knownMSD,1),knownMSD(:,1),'k')
hold on
scatter(1:size(preMSD,1),preMSD(:,1),'g')
hold on
scatter(1:size(MSD,1),MSD(:,1),'b')
hold on
scatter(1:size(QMMSD,1),QMMSD(:,1),'r')
ylim([100,1000])
box on
set(gca,'xscale','log','yscale','log')
xlabel('\tau [frame]');
ylabel('\langledR^{2}\rangle [nm^{2}]');
legend('known','preMSD','4QM','4QM-corrected','Location','northwest')
legend('boxoff')
