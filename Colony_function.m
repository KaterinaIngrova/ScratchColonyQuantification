%% Studium vlivu ropivakain hydrochloridu na migrační a proliferační schopnost nádorových buněk lidského osteosarkomu SaOS-2
% Kateřina Ingrová - BTBIO, 219244
% MATLAB R2020a

results = Colony_function(9.6, 'C:\Users\user\Desktop\diplomka\matlab\colony_EEICT');
%HODNOCENÍ VLIVU LA ROPIVAKAINU NA PROLIFERACI BUNĚK SaOS-2
function [results] = Colony_function(obsah, cesta_k_obrazkum, colony_minimum)

    % 1. Získání seznamu obrázků
    cd(cesta_k_obrazkum);
    images = dir(fullfile(cesta_k_obrazkum, '*.jpg'));

    % 2. Inicializace výsledků
    number_of_images = length(images);
    results = cell(number_of_images, 3);
    results{1,1} = 'Název obrázku';
    results{1,2} = 'Colony Area';
    results{1,3} = 'Počet kolonií';

    % 3. Zpracování jednotlivých obrázků
    for snimek = 1:number_of_images
        % 3.1 Načtení obrázku
        name = images(snimek).name;
        results(snimek+1,1)={name}; %uložení názvu do pole
        image = imread(name);
        figure(1);
        imshow(image);
        disp('------------------------');
        disp(['Název snímku ', name]);
 
        %Zjištění rozměrů obrázku.
        size_of_image = size(image); 
        rows = size_of_image(1);
        columns = size_of_image(2);

        %Vytvoření výřezu snímku kultivační jamky (JE NUTNÉ KLIKNOUT NA MÍSTO
        %KULTIVAČNÍ JAMKY, NASTAVIT ELIPSU A VÝBĚR KULTIVAČNÍ JAMKY POTVRDIT
        %DVOJKLIKEM.

        % Vytvoříme výběr kruhu (interaktivně)
        h_axes = gca;
        kruh = imellipse(h_axes);
        setColor(kruh, 'r'); % Nastavení barvy kruhu na červenou
        wait(kruh); % Čekání na dokončení výběru

        % Získáme pozici kruhu
        position = getPosition(kruh);

        % Uzavřeme figure
        delete(kruh);

        % Získáme střed a poloměr kruhu
        center = [position(1) + position(3)/2, position(2) + position(4)/2];
        radius = min(position(3:4))/2;

        % Vytvoříme masku pro kruh
        [x, y] = meshgrid(1:size(image, 2), 1:size(image, 1));
        maska_kruhu = ((x - center(1)).^2 + (y - center(2)).^2) <= radius^2;

        % Vytvoříme výřez
        vyrez = image;
        for channel = 1:size(image, 3)
            vyrez(:,:,channel) = image(:,:,channel) .* uint8(maska_kruhu);
        end
        
        
        % Získání přiblíženého výřezu bez černých okrajů pomocí bounding boxů
        boundingBox = regionprops(maska_kruhu, 'BoundingBox');
        boundingBox = boundingBox.BoundingBox;

        % Vytvoření finálního vyříznutého snímku
        vyrez_finalni = imcrop(vyrez, boundingBox);

        % Zobrazení finálního vyříznutého snímku
        figure(2);
        subplot(1,3,1)
        imshow(vyrez_finalni);
        title('1');

        % 3.4 %% Zjištění rozměrů vyříznutého snímku a počet černých pixelů reprezentující okraje
        cerny_okraj_pocet_pixelu = sum(vyrez_finalni(:) == 0); %počet černých pixelů výřezu snímků reprezentující okraje
        [rows, cols, ~] = size(vyrez_finalni); 
        celkovy_pocet_pixelu_ve_vyrezu = rows * cols; %celkový počet pixelů výřezu

        % 3.5 %% Zjištění rozměrů k detekci kolonií obsahujících více než 50 buněk
        % najití průměru jamky
        suma_radku = sum(maska_kruhu, 2);
        suma_sloupcu = sum(maska_kruhu,1);
        nejvetsi_suma = max(suma_radku);
        nejvetsi_suma_2 = max(suma_sloupcu);

        % 3,5 cm (35mm)= nejvetsi_suma, 1 mm = 50 buněk
        obsah = obsah;
        prumer_mm = (sqrt(obsah/pi)*2)*10;
        pixel = prumer_mm/nejvetsi_suma_2;
        pocet_pixelu = colony_minimum/(pixel*pixel);
        pocet_pixelu = ceil(pocet_pixelu);

        % Segmentace buněk na základě k-means
        %Převední do formátu double
        img_rgb = vyrez_finalni;
        img_double = double(reshape(img_rgb, [], 3));

        %Určení počtu clusterů
        num_clusters = 3;

        %Aplikace k-means clustering
        [~, centers] = kmeans(img_double, num_clusters);

        %Přiřazení pixelů do clusterů
        [~, cluster_idx] = min(pdist2(img_double, centers), [], 2);

        %Přeformátování clusterového indexu do rozměrů původního obrázku
        segmentedImg = reshape(cluster_idx, size(img_rgb, 1), size(img_rgb, 2));

        %Zjištění počtu pixelů v jednotlivých clusterech
        num_pixels_in_clusters = histcounts(cluster_idx, 1:(num_clusters + 1));

        %Výpis počtu pixelů v jednotlivých clusterech
        %disp('Počet pixelů v jednotlivých clusterech:');
        %disp(num_pixels_in_clusters);

        %Vybrání nejmenšího clusteru - reprezentujicí buňky
        [~, min_cluster_idx] = min(num_pixels_in_clusters);

        %Vytvoření binární masky pro nejmenší cluster
        min_cluster_binary = zeros(size(segmentedImg));
        min_cluster_binary(segmentedImg == min_cluster_idx) = 1;

        %Zobrazení původního a segmentovaného obrázku vedle sebe
        figure(2)
        subplot(1,3,2)
        imshow(min_cluster_binary);
        title('2');
        
        % Aplikace Watershed algoritmu k oddělení jednotlivých kolonií
        %Vytvoření binárního obrazu s použitím negace a otevření (k odstranění
        %malých objektů, šumu a zachování větších objektů)
        binaryImage = ~bwareaopen(~min_cluster_binary, 20);
 
        %Transformace vzdáleností s binárního obrazu
        distanceTransform = -bwdist(~binaryImage);

        %Aplikace algoritmu Watershed
        watershedLabels = watershed(distanceTransform);

        %Jelikož se jedná o přesegmentování, tak následuje úprava
        binaryImage2 = binaryImage;
        binaryImage2(watershedLabels == 0) = 0;

        maskMinima = imextendedmin(distanceTransform,2); %definování extrémních minim

        imposedMinima = imimposemin(distanceTransform,maskMinima); %přidání minim do obrazu D 
        watershedLabels2 = watershed(imposedMinima);
        resultImage = binaryImage;
        resultImage(watershedLabels2 == 0) = 0; %finální obraz po Watershed algoritmu

        
        % 3.7 Zavolání funkce pro spočítání buněk
        [validniBunky, pocetValidnichBunek] = zobrazBunkySCislySVelikosti(resultImage,pocet_pixelu);

        % 3.8 Zobrazení barevného odlišení pravých a nepravých kolonií
         %Zobrazení barevného odlišení pravých a nepravých kolonií

        [rows, columns, ~] = size(validniBunky); %Získání rozměrů obrázku

        novy_obraz = zeros(rows, columns, 3, 'uint8'); % Vytvoření nového obrázku s nulovými hodnotami

        for i = 1:rows
            for j = 1:columns
                % Obarvení buněk na červenou barvu
                if validniBunky(i, j) > 0
                    novy_obraz(i, j, :) = [255, 0, 0];
                end
            end
        end

        binaryImageRGB = im2uint8(cat(3, resultImage, resultImage, resultImage));
        combinedImage = imfuse(binaryImageRGB, novy_obraz,'blend'); %Sjednocení obou snímků pomocí imfuse

        %Extrahování barevných kanálů
        redChannel = combinedImage(:, :, 1);
        greenChannel = combinedImage(:, :, 2);
        blueChannel = combinedImage(:, :, 3);

        %Nalzení šedých pixelů
        blackPixels = redChannel == 128 & greenChannel  == 128 & blueChannel  == 128;
        redPixels = redChannel == 255 & greenChannel  == 128 & blueChannel  == 128;
        pocet_modrych = sum(blackPixels(:) == 1);
        pocet_cervenych = sum(redPixels(:) == 1);

        %Obarvení šedých pixelů na modrou (nepravé kolonie, netvoří 50 buněk)
        redChannel(blackPixels) = 0;
        greenChannel(blackPixels) = 255;
        blueChannel(blackPixels) = 255;

        %Kombinace snímků a získání výsledného odlišovacího snímku
        rgbImage = cat(3, redChannel, greenChannel, blueChannel);
        subplot(1,3,3)
        imshow(rgbImage);
        title('3');

        % 3.9 Výpočet plochy pokryté buňkami
        pomer = pocet_cervenych/(celkovy_pocet_pixelu_ve_vyrezu - cerny_okraj_pocet_pixelu);
        modre_pomer = pocet_modrych/(celkovy_pocet_pixelu_ve_vyrezu - cerny_okraj_pocet_pixelu);
        cerne_pomer = (celkovy_pocet_pixelu_ve_vyrezu - cerny_okraj_pocet_pixelu-pocet_cervenych-pocet_modrych)/(celkovy_pocet_pixelu_ve_vyrezu - cerny_okraj_pocet_pixelu);

        disp(['Poměr pixelu, ktere zaujimaji bunky v jamce: ', num2str(pomer)]);
        disp(['Poměr pixelu, ktere zaujimaji nevalidni bunky v jamce: ', num2str(modre_pomer)]);
        disp(['Poměr pixelu, ktere tvori misku: ', num2str(cerne_pomer)]);
        
        
        % 3.10 Uložení výsledků
        plocha_cervene_cm2 = pomer*obsah;
        plocha_modre_cm2 = modre_pomer*obsah;
        plocha_miska_cm2 = cerne_pomer*obsah;
        disp(['Plocha kolonii: ', num2str(plocha_cervene_cm2), ' cm2']);
        disp(['Plocha nevalidnich kolonii: ', num2str(plocha_modre_cm2), ' cm2']);
        disp(['Plocha misky: ', num2str(plocha_miska_cm2), ' cm2']);
        results(snimek+1,2)={plocha_cervene_cm2};
        results(snimek+1,3)={pocetValidnichBunek};
    end
    
end


% Funkce k počítání kolonií a najití pravých kolonií
function [validniObraz, pocetValidnichBunek]  = zobrazBunkySCislySVelikosti(binarniObraz, minimalniVelikost)

    %Nalezení spojených komponent (bílých oblastí) na binárním obraze
    spojeneKomponenty = bwconncomp(binarniObraz);

    %Získání vlastností kolonií
    vlastnostiBunek = regionprops(spojeneKomponenty, 'Centroid', 'Area', 'PixelList');

    %Filtrace buněk podle minimální velikosti
    validniBunky = vlastnostiBunek([vlastnostiBunek.Area] >= minimalniVelikost);

    %Vytvoření binárního obrazu obsahujícího pouze validní kolonie
    validniObraz = false(size(binarniObraz));

    for i = 1:length(validniBunky)
        pixelList = validniBunky(i).PixelList;
        indices = sub2ind(size(binarniObraz), pixelList(:, 2), pixelList(:, 1));
        validniObraz(indices) = true;
    end
       
    %Počet validních kolonie
    pocetValidnichBunek = length(validniBunky);

    %Zobrazení čísel kolonií vedle středů
    %for i = 1:pocetValidnichBunek
    %    text(validniBunky(i).Centroid(1), validniBunky(i).Centroid(2), num2str(i), 'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold');
    %end

    % Výpis počtu validních kolonií
    fprintf('Počet validních buněk (velikost >= %d): %d\n', minimalniVelikost, pocetValidnichBunek);

end

