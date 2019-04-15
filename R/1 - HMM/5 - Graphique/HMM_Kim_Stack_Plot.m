function [] = HMM_Kim_Stack_Plot(index,mu,sigma,delta,nbRegime_vec,type,n_vec,WorkDataName_vec,smoothProb,Height,Width,heigthDim)
%FONCTION GRAPHIQUE Cette fonction produit le graphique d'une série de
%probabilités lissées
%  Les probabilités lissées sont addictionnées par pas cumulativement de
%  manière à colorer les différentes sections avec la couleur
%  correspondante à l'état
% ———————————————————————————————————————————
% Creation              : 16 janvier 2019
% ———————————————————————————————————————————
% Dernière modification : 6 février  2019
% ———————————————————————————————————————————
colors={'r','y','g','b','c'};

% Longueur du vecteur de la série actuellement affichée
n=n_vec(index);

% Longueur du plus long vecteur parmi les séries qui seront affichées grâce
% à cette fonction
nmax=max(n_vec);

% Nombre de régime de la série actuellement affichée
nbRegime=nbRegime_vec(index);

% Calcul de la largeur nécessaire pour afficher les paramètre du modèle de
% la série actuellement affichée
legendWidth=90+20*(nbRegime-2);
legendHeight=50;

% Initialisation du vecteur pour les probabilités lissées cumulatives
seqKim_stack=zeros(nbRegime+1,n);

% Construction du vecteur par addition cumulative des prob lissées
for i=2:nbRegime+1
    seqKim_stack(i,:)= seqKim_stack(i-1,:)+smoothProb(i-1,:);
end

% Ajout de zéros en début de matrice pour décaller les séries commençant
%	plus tard
if nmax~=n
    seqKim_stack_padded=[zeros(nbRegime+1,nmax-n) seqKim_stack];
else
    seqKim_stack_padded=seqKim_stack;
end

% Initialisation des dimensions et positionnement de la figure en fonction
% du nombre de série à afficher et de l'indice de la série actuellement
% affichée
subplot(Height,Width,[index*4+1,index*4+2,index*4+3]);
tmp=get(gca,'position');
set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])

% Impression des bandes de couleurs correspondantes aux prob lissées
x=[1:nmax, fliplr(1:nmax)];
hold on;
for i=1:nbRegime
	fill(x,...
        [seqKim_stack_padded(i,:),...
        fliplr(seqKim_stack_padded(i+1,:))],...
        string(colors(i)))
end

% Ajout d'une ligne horizontale en y=0 pour faciliter la lecture des
% graphiques
plot(1:nmax,repmat(0.5,1,nmax),'k')

% Modification de l'échelle de l'axe vertical
set(gca,'XLim',[0 nmax],'YLim',[0 1])

% Définition du titre du graphique 
title(strcat(type{index},{' '},'Prob. Lissées P[C_t | X_{1:T}] :',{' '},string(WorkDataName_vec{index})))
hold off;

% Ajout de la légende contenant l'information jugée pertinente pour le
% modèle en fonction de sa structure
h=subplot(Height,Width,index*4+4);
hPos = getpixelposition(h);  

hs = '<html><font size="-2">'; %html start
he = '</font></html>'; %html end

% Si 'HHMM', le vecteur prob contient les probabilités p1-p2-q1-q2
% Si 'HMM',  le vecteur delta contient les probabilités asymptotiques
switch type{index}
    case 'HMM'
        Results={ mu(index,1:nbRegime) ; sigma(index,1:nbRegime); delta(index,1:nbRegime)};
        rowNames={'mu'; 'sigma'; 'delta'};
%         cnh = cellfun(@(x)[hs x he],string(colors(1:nbRegime)),'uni',false); %with html
%         rnh = cellfun(@(x)[hs x he],rowNames,'uni',false); %with html
    case 'HHMM'
        Results={ mu(index,1:nbRegime) ; sigma(index,1:nbRegime); delta(index,1:nbRegime)};
        rowNames={'mu'; 'sigma'; 'delta'};
%         cnh = cellfun(@(x)[hs x he],string(colors(1:nbRegime)),'uni',false); %with html
%         rnh = cellfun(@(x)[hs x he],rowNames,'uni',false); %with html
end
Results=cell2table(Results);
Results.Properties.RowNames = rowNames;

% Insertion des paramètres dans une 'uitable' permettant d'afficher en
% tableau dans un graphique
tabPar = uitable('Data',Results{:,:},...
    'FontSize',8,...
    'ColumnName',string(colors(1:nbRegime)),...
    'RowName',rowNames,...
    'Units', 'pixels',...
    'Position', ([hPos(1) hPos(2)-5 legendWidth+strcmp(string(type{index}),'HHMM')*10 legendHeight]));
set(tabPar, 'ColumnWidth', {legendWidth/(nbRegime-0.5)});
% Les axes sont ensuite effacés pour ne laisser que le tableau
set(h, 'Visible', 'Off')

end

