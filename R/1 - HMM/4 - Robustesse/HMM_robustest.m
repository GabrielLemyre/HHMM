clearvars
% ———————————————————————————————————————————
% MITACS - UdeM - CDPQ
% Cooperation 2018 - 2019
% ———————————————————————————————————————————
% Autheur: Gabriel LEMYRE
% ———————————————————————————————————————————
% Supervision - CDPQ : Océane BENAITEAU
% Supervision - UdeM : Maciej AUGUSTYNIAK
% ———————————————————————————————————————————
% Creation              : 30 novembre 2018
% ———————————————————————————————————————————
% Dernière modification : 14 février  2019
% ———————————————————————————————————————————
% TEST DE ROBUSTESSE
% ———————————————————————————————————————————

% Le dossier suivant doit être synchroniser sur votre ordinateur pour
% pouvoir rouler ce code :
% PVP Direction des risques - udem

% pour le synchroniser, acceder à la page suivante :
% https://cdpcapital.sharepoint.com/sites/Equ-pvprisques/Mesure%20et%20analyse%20quantitative/Forms/AllItems.aspx?RootFolder=%2Fsites%2FEqu%2Dpvprisques%2FMesure%20et%20analyse%20quantitative%2FAnalyse%20%2D%20Global%20Caisse%2FProjets%2Fudem&FolderCTID=0x01200044863C44E8031A4781D4D59A9FE423EB&View=%7BFB93241B%2D1C7F%2D4100%2DA122%2D2FBD0342237F%7D

% Appuyez simplement sur le bouton 'Synchroniser' puis quittez MATLAB et
% ouvrez ce code à nouveau à partir du fichier synchronisé (CDPQ) dans la 
% bar de raccourci de votre desktop.

%% —————————————————————————————————————————
% AJOUT PATH VERS DOSSIERS AUXILIAIRES
%  ————————————————————————————————————————————
breakPath = strsplit(cd,'\');
userPATH = strjoin(breakPath(1:3),'\');

% PATH vers fonctions auxiliaires équipe Risque
addpath(genpath('J:\DEPT\240\APPSDEPT\2210\Gestion du risque\Matlab'));

% PATH vers fonctions auxiliaires HMM - GL
addpath(genpath(strcat(userPATH,'\CDPQ\PVP Direction des risques - udem\Travaux')));

%% —————————————————————————————————————————
% IMPORTATION VARIABLES DU FICHIER DE CONSTRUCTION DES SÉRIES BRUTES
%  ———————————————————————————————————————————
load(strcat(userPATH,'\CDPQ\PVP Direction des risques - udem\Travaux\1 - Base de Données\WorkSpaceBD.mat'))

%% —————————————————————————————————————————
% ENREGISTREMENT DU CHEMIN VERS LE DOSSIER DE FIGURES EXPORTÉES
%  ———————————————————————————————————————————
pathFigures=strcat(userPATH,'\CDPQ\PVP Direction des risques - udem\Travaux\2 - HMM\A - Results\2 - Figures\2 - Robustesse\');

%% —————————————————————————————————————————
%  ////////////////////////////////////////////////////////////////////////
%  SPÉCIFICATION DES PARAMÈTRES
%  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%  ———————————————————————————————————————————
% Quelle aspect doit être testé; produit les séquences d'états et les probabilitées
% 1 - 'startDate'   : Test avec une série mais différentes dates de début
% 2 - 'datasetName' : Test avec séries qui commence toutes à la même date
% 3 - 'nbRegime'    : Test avec une série mais différents nombre de régimes
% 4 - 'typeModele'  : Test avec une série mais différents type de modèles
differenceToTest='typeModele';
%  ———————————————————————————————————————————
% Table contenant toutes les séries
nomTable =       {myTableD myTableM myTableM myTableM myTableM myTableM};

% Si les données n'ont pas toutes la même longueur, la plus longue en
% premier
% datasetName=     {'SPXIndex',   ...
%                   'VIXIndex',  ...
%                   'CPIINDXIndex',...
%                   'CSIBARCIndex',  ...
%                   'MXEFIndex', ...
%                   'USDCADBGNCurncy'};
              
datasetName=     {'SPXIndex',   ...
                  'CCOSTOTIndex',  ...
                  'CPIINDXIndex',...
                  'IPIndex',   ...
                  'CPTICHNGIndex',  ...
                  'CONCCONFIndex'};

% Les dates de départ à comparer doivent être données dans le format
maturite =       [12 ...
                  12 ...
                  12 ...
                  12 ...
                  12 ...
                  12];

% Les dates de départ à comparer doivent être données dans le format
% 'dd-mmm-yyyy'
startDate=       {'01-Jan-2000', ...
                  '01-Jan-2000', ...
                  '01-Jan-2000',...
                  '01-Jan-2000', ...
                  '01-Jan-2000', ...
                  '01-Jan-2000'};

% Sélection de la transformation à apporter à chaque série
transformation = {'logR', ...
                  'logR', ...
                  'logR', ...
                  'logR', ...
                  'logR', ...
                  'logR'};

% Paramètres pour prétraitement des données
freq =           5;

% Pour comparer des séries quotidienne ou hebdomadaire avec mensuel,
% Indiquer 'M' pour toutes les série3s pour garantir qu'elles auront la
% même longueur. Ceci aura pour effet de ne sélectionner que le dernier
% jour ouvrable pour chaque mois
obsFreq =        ['M' 'M' 'M' 'M' 'M' 'M'];

% Indiquer ici si c'est le modèle HMM standard 'HMM' ou le modèle
% hierarchique 'HHMM'.
type =           {'HMM', 'HMM', 'HMM', 'HMM', 'HMM', 'HMM'};
              
% Le nombre de régime dans chaque modèle. 
% Si HHMM, nbRegime = 4 obligatoirement
nbRegime =       [3 3 3 3 3 3];

% Indiquer ici la distribution à utiliser
distribution =   {'Student', 'Student', 'Student', 'Student', 'Student', 'Student'};

%% —————————————————————————————————————————
%  PARAMÈTRES LIÉS À L'AFFICHAGE DES GRAPHIQUES
%  ————————————————————————————————————————————
colors={'r','y','g','b','c'};
nTicksShown=30;
heigthDim=1.1*3/length(nbRegime); % Facteur d'écrasement des graphiques verticalement pour permettre d'afficher les dates

% Spécifie si un graphique de la séries et des séquences doit être produit
% et enregistré lors de l'entraînement du modèle
plotTrue = false; 

%% —————————————————————————————————————————
%  CONSTRUCTION DE L'ÉCHANTILLON à la fréquence demandée et retour de
%                                l'étiquette et du vecteur de dates
%  ———————————————————————————————————————————
% Handle de la fonction prenant l'indice de la série comme entrée
HMM_function=@(x) full_HMM_experience(x,datasetName,distribution,nomTable,maturite,startDate,transformation,freq,nbRegime,plotTrue,pathFigures,obsFreq,type);

% Initialisation des matrices de paramètres
WorkDataName=[];
mu=zeros(length(nbRegime),max(nbRegime));
sigma=zeros(length(nbRegime),max(nbRegime));
delta=zeros(length(nbRegime),max(nbRegime));
prob=zeros(length(nbRegime),4);

% Entrainement des modèles
[data1,dates1,mu(1,1:nbRegime(1)),sigma(1,1:nbRegime(1)),delta(1,1:nbRegime(1)),~,~,sequence_Lisse1,smoothProb1,p_ct_x1tm1_1,WorkDataName{1}]= HMM_function(1);
[data2,dates2,mu(2,1:nbRegime(2)),sigma(2,1:nbRegime(2)),delta(2,1:nbRegime(2)),~,~,sequence_Lisse2,smoothProb2,p_ct_x1tm1_2,WorkDataName{2}]= HMM_function(2);
[data3,dates3,mu(3,1:nbRegime(3)),sigma(3,1:nbRegime(3)),delta(3,1:nbRegime(3)),~,~,sequence_Lisse3,smoothProb3,p_ct_x1tm1_3,WorkDataName{3}]= HMM_function(3);
[data4,dates4,mu(4,1:nbRegime(4)),sigma(4,1:nbRegime(4)),delta(4,1:nbRegime(4)),~,~,sequence_Lisse4,smoothProb4,p_ct_x1tm1_4,WorkDataName{4}]= HMM_function(4);
[data5,dates5,mu(5,1:nbRegime(5)),sigma(5,1:nbRegime(5)),delta(5,1:nbRegime(5)),~,~,sequence_Lisse5,smoothProb5,p_ct_x1tm1_5,WorkDataName{5}]= HMM_function(5);
[data6,dates6,mu(6,1:nbRegime(6)),sigma(6,1:nbRegime(6)),delta(6,1:nbRegime(6)),~,~,sequence_Lisse6,smoothProb6,p_ct_x1tm1_6,WorkDataName{6}]= HMM_function(6);

% Obtention des longueurs des diverses séries
n_vec=[length(data1) ...
       length(data2) ... 
       length(data3) ...
       length(data4) ... 
       length(data5) ... 
       length(data6)];

[nmax,indmax]=max(n_vec);

% Changement format des dates
formatOut = 'yy';
dates=datestr(dates1,formatOut);


%% —————————————————————————————————————————
% CONSTRUCTION GRAPHIQUES DES PROBABILITÉS LISSÉES
%  ———————————————————————————————————————————
% Hauteur et largeur de la figure en nombre de sections
h_sub=length(datasetName)+1;
w_sub=4;

% Distance entre points à mettre sur l'axe horizontale ainsi que les
% emplacements de ceux-ci
nTicks=floor(nmax/nTicksShown);
xTicksValues=1:nTicks:nmax;

%  ——————————————————————————
%  ——————————————————————————

% Graphique de la série prétraitée
subplot(h_sub,w_sub,[1,2,3]);
tmp=get(gca,'position');
set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
hold on;
plot(1:length(dates1),data1)
set(gca,'XLim',[0 nmax],'YLim',[-1.2*max(abs(data1)) 1.2*max(abs(data1))])
title(string(datasetName{1}))
% Ajout d'une ligne horizontale en y=0 pour faciliter la lecture
plot(1:nmax,zeros(1,nmax),'k')
hold off;
h=subplot(h_sub,w_sub,4);
set(h, 'Visible', 'Off')

%  ——————————————————————————
%  ——————————————————————————
Plot_func = @(x,y) HMM_Kim_Stack_Plot(x,mu,sigma,delta,nbRegime,type,n_vec,WorkDataName,y,h_sub,w_sub,heigthDim);

Plot_func(1,smoothProb1)
Plot_func(2,smoothProb2)
Plot_func(3,smoothProb3)
Plot_func(4,smoothProb4)
Plot_func(5,smoothProb5)
Plot_func(6,smoothProb6)
%  ——————————————————————————
%  ——————————————————————————

set(findobj(gcf,'type','axes'),'FontName','Arial',...
    'FontSize',8,...
    'LineWidth', 0.5,...
    'TickLength',[0.005, 0.005],...
    'XTick', xTicksValues,...
    'XTickLabel', dates(1:nTicks:nmax,:),...
    'XTickLabelRotation',90);

% 
set(gcf,'PaperUnits', 'inches', ...
    'PaperPosition',[0, 0, 20, 10], ...
    'Units', 'pixels',...
    'papersize', [1800 1200])

switch differenceToTest
    case 'startDate'
        s=strcat(startDate,'-');
    case 'datasetName'
        s=strcat(datasetName,'-');
    case 'nbRegime'
        s=strcat(string(nbRegime),'-');
    case 'typeModele'
        s=strcat('(',string(type),'-',string(nbRegime),')');
end

print(strcat(pathFigures,...
        WorkDataName{1},'_',...
        string(nbRegime(1)),'_regimes_Test_Robustesse_',...
        strcat(s{:}),...
        '_KIM.png'),...
    '-dpng','-r500')

close(gcf)


%% —————————————————————————————————————————
% CONSTRUCTION GRAPHIQUES DES PROBABILITÉS FILTRÉES
%  ———————————————————————————————————————————   
% if n1~=n2
%     sequence2_padded=[1 zeros(1,n1-n2-1)*NaN sequence2];
% else
%     sequence2_padded=sequence2;
% end
% if n1~=n3
%     sequence3_padded=[1 zeros(1,n1-n3-1)*NaN sequence3];
% else
%     sequence3_padded=sequence3;
% end
% if n1~=n4
%     sequence4_padded=[1 zeros(1,n1-n4-1)*NaN sequence4];
% else
%     sequence4_padded=sequence4;
% end
% 
% subplot(5,4,[1,2,3]);
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
% plot(data1)
% ylim([-1.2*max(abs(data1)) 1.2*max(abs(data1))])
% title(strcat(datasetName{1},' with ',{' '},string(nbRegime(1)),' Regimes'))
% 
% subplot(5,4,[5,6,7]);
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
% plot(1:n1,sequence1)
% ylim([0.95 (nbRegime(1)+0.05)])
% % Définition du titre du graphique VITERBI 1
% switch differenceToTest
% 	case 'startDate'
%         title(strcat('Séquence Viterbi - date début :',{' '},string(startDate{1})))
%     case 'datasetName'
%         title(strcat('Séquence Viterbi - Séries de :',{' '},string(datasetName{1})))
%     case 'nbRegime'
%         title(strcat('Séquence Viterbi - nbRegime :',{' '},string(nbRegime(1))))
% end
% 
% h=subplot(5,4,8);
% hPos = get(h, 'Position');  
% Results={ mu1 ; sigma1};
% rowNames={'mu_1'; 'sigma_1'};
% Results=cell2table(Results);
% Results.Properties.RowNames = rowNames;
% 
% uitable('Data',Results{:,:},...
%     'ColumnName',string(colors(1:nbRegime(1))),...
%     'RowName',Results.Properties.RowNames,...
%     'Units', 'Normalized',...
%     'Position', ([hPos(1:2) legendWidth_1 legendHeight]));
% 
% set(h, 'Visible', 'Off')
% 
% subplot(5,4,[9,10,11]);
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
% plot(1:n1,sequence2_padded)
% ylim([0.95 (nbRegime(2)+0.05)])
% % Définition du titre du graphique VITERBI 2
% switch differenceToTest
% 	case 'startDate'
%         title(strcat('Séquence Viterbi - date début :',{' '},string(startDate{2})))
%     case 'datasetName'
%         title(strcat('Séquence Viterbi - Séries de :',{' '},string(datasetName{2})))
%     case 'nbRegime'
%         title(strcat('Séquence Viterbi - nbRegime :',{' '},string(nbRegime(2))))
% end
% 
% h=subplot(5,4,12);
% hPos = get(h, 'Position');  
% Results={ mu2 ; sigma2};
% rowNames={'mu_2'; 'sigma_2'};
% Results=cell2table(Results);
% Results.Properties.RowNames = rowNames;
% 
% uitable('Data',Results{:,:},...
%     'ColumnName',string(colors(1:nbRegime(2))),...
%     'RowName',Results.Properties.RowNames,...
%     'Units', 'Normalized',...
%     'Position', ([hPos(1:2) legendWidth_2 legendHeight]));
% 
% set(h, 'Visible', 'Off')
% 
% subplot(5,4,[13,14,15]);
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
% plot(1:n1,sequence3_padded)
% ylim([0.95 (nbRegime(3)+0.05)])
% % Définition du titre du graphique VITERBI 3
% switch differenceToTest
% 	case 'startDate'
%         title(strcat('Séquence Viterbi - date début :',{' '},string(startDate{3})))
%     case 'datasetName'
%         title(strcat('Séquence Viterbi - Séries de :',{' '},string(datasetName{3})))
%     case 'nbRegime'
%         title(strcat('Séquence Viterbi - nbRegime :',{' '},string(nbRegime(3))))
% end
% 
% h=subplot(5,4,16);
% hPos = get(h, 'Position');  
% Results={ mu3 ; sigma3};
% rowNames={'mu_3'; 'sigma_3'};
% Results=cell2table(Results);
% Results.Properties.RowNames = rowNames;
% 
% uitable('Data',Results{:,:},...
%     'ColumnName',string(colors(1:nbRegime(3))),...
%     'RowName',Results.Properties.RowNames,...
%     'Units', 'Normalized',...
%     'Position', ([hPos(1:2) legendWidth_3 legendHeight]));
% 
% set(h, 'Visible', 'Off')
% 
% subplot(5,4,[17,18,19]);
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
% plot(1:n1,sequence4_padded)
% ylim([0.95 (nbRegime(4)+0.05)])
% % Définition du titre du graphique VITERBI 4
% switch differenceToTest
% 	case 'startDate'
%         title(strcat('Séquence Viterbi - date début :',{' '},string(startDate{4})))
%     case 'datasetName'
%         title(strcat('Séquence Viterbi - Séries de :',{' '},string(datasetName{4})))
%     case 'nbRegime'
%         title(strcat('Séquence Viterbi - nbRegime :',{' '},string(nbRegime(4))))
% end
% 
% h=subplot(5,4,20);
% hPos = get(h, 'Position');  
% Results={ mu4 ; sigma4};
% rowNames={'mu_4'; 'sigma_4'};
% Results=cell2table(Results);
% Results.Properties.RowNames = rowNames;
% 
% uitable('Data',Results{:,:},...
%     'ColumnName',string(colors(1:nbRegime(4))),...
%     'RowName',Results.Properties.RowNames,...
%     'Units', 'Normalized',...
%     'Position', ([hPos(1:2) legendWidth_4 legendHeight]));
% 
% set(h, 'Visible', 'Off')
% 
% set(findobj(gcf,'type','axes'),'FontName','Arial',...
%     'FontSize',9,...
%     'LineWidth', 0.5,...
%     'TickLength',[0.005, 0.005],...
%     'XTick', xTicksValues,...
%     'XTickLabel', dates(1:nTicks:n1,:),...
%     'XTickLabelRotation',90);
% 
% set(gcf,'PaperPositionMode','auto',...
%     'Units', 'Normalized',...
%     'OuterPosition', [0, 0.04, 1, 0.96])
% 
% s=strcat(etiquette_operation,'-');
% 
% 
% print(strcat(pathFigures,...
%         WorkDataName1,'_',...
%         string(nbRegime(1)),'_regimes_Test_Robustesse_',...
%         strcat(s{:}),...
%         '_VITERBI.png'),...
%     '-dpng','-r500')
% 
% close(gcf)












