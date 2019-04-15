function [] = HMM_Kim_Stack_Plot(index,mu,sigma,delta,nbRegime_vec,type,n_vec,WorkDataName_vec,smoothProb,Height,Width,heigthDim)
%FONCTION GRAPHIQUE Cette fonction produit le graphique d'une s�rie de
%probabilit�s liss�es
%  Les probabilit�s liss�es sont addictionn�es par pas cumulativement de
%  mani�re � colorer les diff�rentes sections avec la couleur
%  correspondante � l'�tat
% �������������������������������������������
% Creation              : 16 janvier 2019
% �������������������������������������������
% Derni�re modification : 6 f�vrier  2019
% �������������������������������������������
colors={'r','y','g','b','c'};

% Longueur du vecteur de la s�rie actuellement affich�e
n=n_vec(index);

% Longueur du plus long vecteur parmi les s�ries qui seront affich�es gr�ce
% � cette fonction
nmax=max(n_vec);

% Nombre de r�gime de la s�rie actuellement affich�e
nbRegime=nbRegime_vec(index);

% Calcul de la largeur n�cessaire pour afficher les param�tre du mod�le de
% la s�rie actuellement affich�e
legendWidth=90+20*(nbRegime-2);
legendHeight=50;

% Initialisation du vecteur pour les probabilit�s liss�es cumulatives
seqKim_stack=zeros(nbRegime+1,n);

% Construction du vecteur par addition cumulative des prob liss�es
for i=2:nbRegime+1
    seqKim_stack(i,:)= seqKim_stack(i-1,:)+smoothProb(i-1,:);
end

% Ajout de z�ros en d�but de matrice pour d�caller les s�ries commen�ant
%	plus tard
if nmax~=n
    seqKim_stack_padded=[zeros(nbRegime+1,nmax-n) seqKim_stack];
else
    seqKim_stack_padded=seqKim_stack;
end

% Initialisation des dimensions et positionnement de la figure en fonction
% du nombre de s�rie � afficher et de l'indice de la s�rie actuellement
% affich�e
subplot(Height,Width,[index*4+1,index*4+2,index*4+3]);
tmp=get(gca,'position');
set(gca,'position',[tmp(1) tmp(2) tmp(3) heigthDim*tmp(4)])

% Impression des bandes de couleurs correspondantes aux prob liss�es
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

% Modification de l'�chelle de l'axe vertical
set(gca,'XLim',[0 nmax],'YLim',[0 1])

% D�finition du titre du graphique 
title(strcat(type{index},{' '},'Prob. Liss�es P[C_t | X_{1:T}] :',{' '},string(WorkDataName_vec{index})))
hold off;

% Ajout de la l�gende contenant l'information jug�e pertinente pour le
% mod�le en fonction de sa structure
h=subplot(Height,Width,index*4+4);
hPos = getpixelposition(h);  

hs = '<html><font size="-2">'; %html start
he = '</font></html>'; %html end

% Si 'HHMM', le vecteur prob contient les probabilit�s p1-p2-q1-q2
% Si 'HMM',  le vecteur delta contient les probabilit�s asymptotiques
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

% Insertion des param�tres dans une 'uitable' permettant d'afficher en
% tableau dans un graphique
tabPar = uitable('Data',Results{:,:},...
    'FontSize',8,...
    'ColumnName',string(colors(1:nbRegime)),...
    'RowName',rowNames,...
    'Units', 'pixels',...
    'Position', ([hPos(1) hPos(2)-5 legendWidth+strcmp(string(type{index}),'HHMM')*10 legendHeight]));
set(tabPar, 'ColumnWidth', {legendWidth/(nbRegime-0.5)});
% Les axes sont ensuite effac�s pour ne laisser que le tableau
set(h, 'Visible', 'Off')

end

