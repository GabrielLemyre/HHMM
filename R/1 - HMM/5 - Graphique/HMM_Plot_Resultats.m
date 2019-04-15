function [] = HMM_Plot_Resultats(TypeModele,barGraphTrue, OrientedGraph ,datasetName, Data,dates,mu, sigma,Distribution,gamma,delta,prob,tempsMoyenSejour,seqKim, pathFigures)
%FONCTION GRAPHIQUE Fonction permettant de produire les graphiques
%   Cette fonction produit les graphiques de la séquence d'observations, de
%   la séquence d'état la plus probable avec algorithme de Viterbi ainsi
%   que les probabilitées lissées avec filtre de Kim. Les paramètres du
%   modèles final sont aussi imprimmé au bas de la figure
% ———————————————————————————————————————————
% Paramètres :
% 	      TypeModele,       : 
% 	      barGraphTrue,     :
% 	      OrientedGraph,    :
% 	      datasetName,      :
% 	      Data,             :
% 	      dates,            :
% 	      mu,               :
% 	      sigma,            :
% 	      Distribution,     :
% 	      gamma,            :    
% 	      delta,            :
% 	      prob,             :
% 	      tempsMoyenSejour, :
% 	      seqKim,           :
% 	      pathFigures       :
% ———————————————————————————————————————————
% Creation              : 18 décembre 2018
% ———————————————————————————————————————————
% Dernière modification : 25 février 2019
% ———————————————————————————————————————————


% ———————————————————————————————————————————
% PARAMÈTRES SUR LES DONNÉES ET LES MODÈLES
% ———————————————————————————————————————————
% Longueur du vecteur de données
n=length(Data);

% Nombre de régimes dans le modèle
nbRegime=size(seqKim,1);

% Choix du nombre de degré de liberté pour la distribution de Student
% Aussi défini dans : normal_HMM_HamiltonFilter -- normal_HMM_mle -- normal_HMM_Viterbi
nu=4;


% ———————————————————————————————————————————
% PARAMÈTRES POUR L'AFFICHAGE ET LES GRAPHIQUES
% ———————————————————————————————————————————

% Définition du nombre de colonnes et de lignes du quadrillage d'affichage
% des graphiques
subRow = 5;
subColumn = 4;

% Dimension des boîtes contenant les graphiques
legendWidth=85+25*(nbRegime-2);
legendHeight=45;
heigthDim=0.72; % Facteur d'écrasement des graphiques verticalement pour permettre d'afficher les dates

% Initialisation de l'indice correspondant au numéro de position des
% graphiques (voir appel de nouveau subplots)
index_plot = 0;

% Largeur du trait dans les graphiques de distributions statistiques
plotLineWidth = 1.5;

% Nombre de points à inscrire sur l'axe horizontale
nTicksShown=min(n,30);

% Nombre d'observations entre chaque point à afficher sur l'axe horizontale
nTicks=floor(n/nTicksShown);

% Construction du vecteur contenant les indices à afficher sur l'axe
% horizontale
xTicksValues=1:nTicks:n;

% Puisque le code pour les graphiques orientés ne fonctionnent qu'avec 4
% régimes pour l'instant, ceci affiche un message d'erreur dans le cas où
% le nombre de régimes n'est pas 4 et qu'on demande l'affichage du
% graphique orienté dans l'appel de la fonction
if OrientedGraph
    if nbRegime ~= 4
        error(sprintf("Erreur. \nLes graphiques orientées du processus de changement de régimes ne peuvent être imprimées que si le modèle comporte 4 régimes. \nLe graphique sera remplacé par l'affichage habituel. GL"));
        OrientedGraph=false;
    end
end

% Modification des dimensions pour l'affichage de résultats sur HHMM
if (strcmp(TypeModele,'HHMM')) && (~OrientedGraph)
    % Si HHMM, ajout d'une ligne pour afficher les probabilités lissées pour
    % les régimes parents
    subRow = subRow + 1;
	% Réduction hauteur graphique pour laisser place aux parents
    heigthDim = heigthDim*0.9;
end

% Choix de la couleur de l'affichage en fonction du nombre de régimes
switch nbRegime
    case 2
        colors_Names={'Rouge','Bleu'};
        colors=[[1    0    0]; ...
                [0,         0,    0.6000]];
    case 3
        colors_Names={'Rouge','Cyan','Bleu'};
        colors=[[1    0    0]; ...
                [0.2000    0.6000    1.0000]; ...
                [0,         0,    0.6000]];
    case 4
        colors_Names={'Rouge','Saumon','Cyan','Bleu'};
        colors=[[1    0    0]; ...
                [0.9804,    0.5020,    0.4471]; ...
                [0.2000    0.6000    1.0000]; ...
                [0,         0,    0.6000]];
    otherwise
        colors_Names={'Rouge','Orange','Vert','Bleu','Noir'};
        colors=[[1    0    0]; ...
                [0.9290, 0.6940, 0.1250]; ...
                [0,1,0]; ...
                [0, 0, 1]; ...
                [0.25, 0.25, 0.25]];
end




% ———————————————————————————————————————————
% CONSTRUCTION VECTEUR DE PROBABILITÉS ET DE SÉQUENCE ÉTATS PLUS PROBABLE
% ———————————————————————————————————————————
% Construction de la séquence d'états la plus probable a posteriori
sequence_plus_probable_etats = zeros(1,n);
for i=n:-1:1
    [~,whichMax]=max(seqKim(:,i));
    sequence_plus_probable_etats(i)=whichMax;
end

numberReg = string([sum(sequence_plus_probable_etats==1), sum(sequence_plus_probable_etats==2), sum(sequence_plus_probable_etats==3), sum(sequence_plus_probable_etats==4)]);
numberRegNUM = str2double(numberReg);
nbRegimeRealiser = sum(numberRegNUM>0);

% Somme cumulative des probabilités lissées
seqKim_stack=zeros(nbRegime+1,n);
for i=2:nbRegime+1
	seqKim_stack(i,:)= seqKim_stack(i-1,:)+seqKim(i-1,:);
end

%  ————————————————————————————————————————————
%  GRAPHIQUE DE LA SÉRIE TELLE QU'UTILISÉE POUR L'ENTRAINEMENT
%  ————————————————————————————————————————————
subplot(subRow,subColumn,(index_plot + (1:4)));

if barGraphTrue
    bar(Data, 'FaceColor','flat','EdgeColor','flat','CData',colors(sequence_plus_probable_etats,:));
    set(gca,'XLim',[0 n],'YLim',[-1.2*max(abs(Data)) 1.2*max(abs(Data))])
else
    plot(Data)
    set(gca,'XLim',[0 n],'YLim',[-1.2*max(abs(Data)) 1.2*max(abs(Data))])
end

title(strcat(TypeModele,' -- ',{' '},datasetName,' with ',{' '},string(nbRegime),' Regimes'))
tmpBarPlot=get(gca,'position');
set(gca,'position',[tmpBarPlot(1) tmpBarPlot(2) tmpBarPlot(3) heigthDim*tmpBarPlot(4)])


if OrientedGraph
    set(findobj(gcf,'type','axes'),'FontName','Arial',...
        'FontSize',9,...
        'LineWidth', 0.5,...
        'TickLength',[0.005, 0.005],...
        'XTick', xTicksValues,...
        'XTickLabel', dates(1:nTicks:n,:),...
        'XTickLabelRotation',90);
end

index_plot = index_plot + 4;
%  ————————————————————————————————————————————
%  GRAPHIQUE SÉQUENCE D'ÉTAT LA PLUS PROBABLE (VITERBI)
%  ————————————————————————————————————————————
% subplot(4,2,[3 4]);
% plot(seqViterbi)
% ylim([0.95 (nbRegime+0.05)])
% title('Sequence Viterbi')
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])

% Si l'appel de la fonction demande qu'un graphique orienté du modèle soit
% produit, celui-ci remplace les probabilitées lissées ainsi que les
% distribution et boxplots
if OrientedGraph
    m=subplot(subRow,subColumn,(index_plot + (1:12)));
    mPos = getpixelposition(m);   
    if strcmp(TypeModele,'HHMM')
        HMM_orientedGraph(TypeModele,gamma,prob,datasetName,pathFigures,1,false,true);
        title('Graphique Hiérarchie')
    else
        HMM_orientedGraph(TypeModele,gamma,prob,datasetName,pathFigures,1,false,false);
        title('Processus Changement de Régime')
    end
    tmp=get(gca,'position');
    set(gca,'position',[tmp(1) tmp(2)-0.045 tmp(3) tmp(4)])
    % set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])
    
    index_plot = index_plot + 12;
else
    %  ————————————————————————————————————————————
    %  GRAPHIQUE DE LA SÉQUENCE DE PROBABILITÉS LISSÉES (FILTRE DE KIM)
    %  ————————————————————————————————————————————
    m=subplot(subRow,subColumn,(index_plot + (1:4)));
    mPos = getpixelposition(m);   
    x=[1:n, fliplr(1:n)];
    hold on;
    for i=1:nbRegime
        fill(x,...
            [seqKim_stack(i,:),...
            fliplr(seqKim_stack(i+1,:))],...
            colors(i,:))
    end
    plot(1:n,repmat(0.5,1,n),'k')
    set(gca,'XLim',[0 n],'YLim',[0 1])
    title('Séquence Probabilités lissées P[C_t | X_1:T]')
    hold off;
    tmp=get(gca,'position');
    set(gca,'position',[tmp(1) tmp(2)-0.015*strcmp(TypeModele,'HHMM') tmpBarPlot(3) heigthDim*tmp(4)])
    % set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])

    set(findobj(gcf,'type','axes'),'FontName','Arial',...
        'FontSize',9,...
        'LineWidth', 0.5,...
        'TickLength',[0.005, 0.005],...
        'XTick', xTicksValues,...
        'XTickLabel', dates(1:nTicks:n,:),...
        'XTickLabelRotation',90);

    index_plot = index_plot + 4;


    %  ————————————————————————————————————————————
    %  GRAPHIQUE DE LA SÉQUENCE DE PROBABILITÉS LISSÉES (FILTRE DE KIM) DES
    %  RÉGIMES PARENTS
    %  ————————————————————————————————————————————
    if strcmp(TypeModele,'HHMM')
        m=subplot(subRow,subColumn,(index_plot + (1:4)));
        mPos = getpixelposition(m); 

        % Addition des probabilités lissées par régime parent
        seqKim_parent=[seqKim(1,:)+seqKim(2,:) ; seqKim(3,:)+seqKim(4,:)];

        % Somme cumulative des probabilités lissées par régime parent
        seqKim_stack_parent=zeros(3,n);
        for i=2:3
            seqKim_stack_parent(i,:)= seqKim_stack_parent(i-1,:)+seqKim_parent(i-1,:);
        end

        % Affichage coloré des probabilités lissées cumulatives par parent
        colors_HHMM = [[0,0,0] ; ...
                       [1,1,1]];
        x=[1:n, fliplr(1:n)];
        hold on;
        for i=1:2
            fill(x,...
                [seqKim_stack_parent(i,:),...
                fliplr(seqKim_stack_parent(i+1,:))],...
                colors_HHMM(i,:))
        end
        plot(1:n,repmat(0.5,1,n),'b')
        set(gca,'XLim',[0 n],'YLim',[0 1])
        title('Séquence Probabilités lissées P[C_t | X_1:T] : Régimes Parents')
        hold off;
        tmp=get(gca,'position');
        set(gca,'position',[tmp(1) tmp(2)-0.03 tmpBarPlot(3) heigthDim*tmp(4)])

        set(findobj(gcf,'type','axes'),'FontName','Arial',...
            'FontSize',9,...
            'LineWidth', 0.5,...
            'TickLength',[0.005, 0.005],...
            'XTick', xTicksValues,...
            'XTickLabel', dates(1:nTicks:n,:),...
            'XTickLabelRotation',90);

        index_plot = index_plot + 4;
    end


    %  ————————————————————————————————————————————
    %  GRAPHIQUE DISTRIBUTIONS OBTENUES
    %  ————————————————————————————————————————————
    subplot(subRow,subColumn,(index_plot + [1,2,5,6]));
    muRealiser = mu(numberRegNUM>0);
    sigmaRealiser = sigma(numberRegNUM>0);
    pas = 2*(max(muRealiser)+3*max(sigmaRealiser))/100000;
    colorsRealiser = colors(numberRegNUM>0,:);
    switch Distribution
        case 'Student'
            x = -(max(muRealiser)+3*max(sigmaRealiser)):pas:(max(muRealiser)+3*max(sigmaRealiser));
            y = zeros(nbRegimeRealiser,length(x));
            for i=1:nbRegimeRealiser
                y(i,:) = tpdf((x-muRealiser(i))/sigmaRealiser(i),nu)./sigmaRealiser(i);
            end
            plot(x,y(1,:),'color',colorsRealiser(1,:),'LineStyle','-','LineWidth',plotLineWidth)
            hold on
            for i=2:nbRegimeRealiser
            	plot(x,y(i,:),'color',colorsRealiser(i,:),'LineStyle','-','LineWidth',plotLineWidth)
            end
            legend(strcat('mu = ',{' '},string(round(muRealiser,3)),{' '},'sigma = ',{' '},string(round(sigmaRealiser,3))).','Location','best')
            hold off

        case 'Normale'
            x = -(max(muRealiser)+3*max(sigmaRealiser)):pas:(max(muRealiser)+3*max(sigmaRealiser));
            y = zeros(nbRegimeRealiser,length(x));
            for i=1:nbRegimeRealiser
                y(i,:) = normpdf(x,muRealiser(i),sigmaRealiser(i));
            end
            plot(x,y(1,:),'color',colorsRealiser(1,:),'LineStyle','-','LineWidth',plotLineWidth)
            hold on
            for i=2:nbRegimeRealiser
            	plot(x,y(i,:),'color',colorsRealiser(i,:),'LineStyle','-','LineWidth',plotLineWidth)
            end
            legend(strcat('mu = ',{' '},string(round(muRealiser,3)),{' '},'sigma = ',{' '},string(round(sigmaRealiser,3))).','Location','best')
            hold off
    end

    tmp=get(gca,'position');
    set(gca,'position',[tmp(1)+0.005 tmp(2) tmp(3) heigthDim*tmp(4)])


    %  ————————————————————————————————————————————
    %  GRAPHIQUES BOX PLOT PAR RÉGIME
    %  ————————————————————————————————————————————
    subplot(subRow,subColumn,(index_plot + [3,4,7,8]));
    boxplot(Data,sequence_plus_probable_etats)
    xlabel('Régimes-(fréquence)');
    ylabel('Rendements');
    indx = str2num(cell2mat(xticklabels));
    xticklabels(strcat(colors_Names(indx),'-(',string(round(100*(str2double(numberReg(indx))/sum(str2double(numberReg))),2)),'%)'));
    title('Boxplot par régime');

    tmp1=get(gca,'position');
    set(gca,'position',[tmp1(1) tmp(2) tmp(3)*0.8 heigthDim*tmp(4)])

    index_plot = index_plot + 8;
end
%  ————————————————————————————————————————————
%  TABLEAU DES PARAMÈTRES DU MODÈLE
%  ————————————————————————————————————————————

h=subplot(subRow,subColumn,(index_plot + [1,2]));
hPos = getpixelposition(h);   
Results=[ string(mu) ; string(sigma) ; string(delta); string(tempsMoyenSejour.')];

%  ——————————————————————————————
    % Fonction permettant de colorer les cellules associée à la
    % plus haute/basse volatilité et haute/basse moyenne
    colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
    % Extraction de l'indice de la colonne d'intérêt
    [~,WhichVolMax] = max(sigma);
    [~,WhichVolMin] = max(-sigma);
    [~,WhichMeanMax] = max(mu);
    [~,WhichMeanMin] = max(-mu);

    % Ajout code HTML permettant de colorer en 'jaune' les cellules une à une
    volMax = colergen('#FFFF00',Results(2,WhichVolMax));
    Results(2,WhichVolMax)=strcat(volMax{:});
    
    volMin = colergen('#32CD32',Results(2,WhichVolMin));
    Results(2,WhichVolMin)=strcat(volMin{:});
    
    meanMax = colergen('#32CD32',Results(1,WhichMeanMax));
    Results(1,WhichMeanMax)=strcat(meanMax{:});
    
    meanMin = colergen('#FFFF00',Results(1,WhichMeanMin));
    Results(1,WhichMeanMin)=strcat(meanMin{:});
    
%  ——————————————————————————————

Results=cellstr(num2cell(Results));
rowNames={'mu'; 'sigma'; 'delta'; 'TMS'};

uitable('Data',Results,...
    'ColumnName',string(colors_Names(1:nbRegime)),...
    'RowName',rowNames,...
    'Position', ([hPos(1) 20 legendWidth legendHeight]));

set(h, 'Visible', 'Off')

%  ————————————————————————————————————————————
%  TABLEAU DE LA MATRICE DE TRANSITION DU MODÈLE
%  ————————————————————————————————————————————
h=subplot(subRow,subColumn,(index_plot + [3,4]));
hPos = getpixelposition(h);  
switch TypeModele
    case 'HMM'
        RowNames=cellstr(strcat('etat',string(1:nbRegime)));
    case 'HHMM'
        RowNames=['S_{1}'; 'S_{2}'; 'S_{3}'; 'S_{4}' ];
end


Results=cellstr(strcat(string(round(100*gamma,2)),'%'));

uitable('Data',Results,...
    'ColumnName',string(colors_Names(1:nbRegime)),...
    'RowName',RowNames,...
    'Units', 'pixels',...
    'Position', ([hPos(1)+20 40-10*(nbRegime-2) legendWidth legendHeight+10*(nbRegime-4)]));


set(h, 'Visible', 'Off')

%  ————————————————————————————————————————————
%  PARAMÈTRES LIÉS À L'AFFICHAGE DES GRAPHIQUES
%  ————————————————————————————————————————————

set(gcf,'PaperUnits', 'inches', ...
    'PaperPosition',[0, 0, 17, 10], ...
    'Units', 'pixels',...
    'papersize', [1400 900])


% Creation du nom du fichier a créer pour recevoir les résultats
newSubFolder = char(strcat(pathFigures,datasetName));

% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

if OrientedGraph
    switch TypeModele
        case 'HHMM'
            print(strcat(newSubFolder,'\',datasetName,'_Modele_',TypeModele,'_Regime_Switching_Graph.png'),'-dpng','-r500')
        case 'HMM'
            print(strcat(newSubFolder,'\',datasetName,'_Modele_',string(nbRegime),'_regimes_',TypeModele,'_Regime_Switching_Graph.png'),'-dpng','-r500')
    end
else
    switch TypeModele
        case 'HHMM'
            print(strcat(newSubFolder,'\',datasetName,'_Modele_',TypeModele,'.png'),'-dpng','-r500')
        case 'HMM'
            print(strcat(newSubFolder,'\',datasetName,'_Modele_',string(nbRegime),'_regimes_',TypeModele,'.png'),'-dpng','-r500')
    end
end

close(gcf)
end

