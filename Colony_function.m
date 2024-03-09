%Effect of gentamicin on colony formation in the SaOS-2 cell line

%Example of calling function with the obtained results
%results = Colony_function(9.6, 'C:\Users\user\Desktop\diplomka\matlab\colony_20_11_nahore', 0.03);

%Evaluation function - inputs are the path to the images, the area of the well in cm2 and the minimum colony size. 
%The obtained results are stored in the cell (number of colonies, area covered by colonies in cm2).
function [results] = Colony_function(well_area, path_to_image, colony_minimum)

    %List of images 
    cd(path_to_image);
    files = dir(fullfile(path_to_image, '*.*'));
    images = files(~[files.isdir]);

    %Inicialization of cell for storing results
    number_of_images = length(images);
    results = cell(number_of_images, 3);
    results{1,1} = 'Image name';
    results{1,2} = 'Colony area';
    results{1,3} = 'Colony count';

    %For loop for loading current image
    for i = 1:number_of_images
        name = images(i).name; %name of image
        results(i+1,1)={name}; %saving the name of image to created cell
        image = imread(name); %load concrete image
        
        %Show the image for individual marking well
        figure(1);
        imshow(image);
        disp('------------------------');
        disp(['Název snímku ', name]);
 
        %Determination of image size
        size_of_image = size(image); 
        rows = size_of_image(1);
        columns = size_of_image(2);

        %Create a slice of the image of the cultivation well (MUST CLICK ON THE PLACE
        %CULTURE WELL, SET THE ELLIPSE AND CONFIRM THE SELECTION OF THE
        %CULTURE WELL DOUBLE CLICK)
        h_axes = gca;
        circle = imellipse(h_axes);
        setColor(circle, 'r');
        wait(circle);

        %Get position of circle
        position = getPosition(circle);

        %Close figure
        delete(circle);

        %Get the centre and radius of the circle
        center = [position(1) + position(3)/2, position(2) + position(4)/2];
        radius = min(position(3:4))/2;

        %Create mask for circle
        [x, y] = meshgrid(1:size(image, 2), 1:size(image, 1));
        mask_circle = ((x - center(1)).^2 + (y - center(2)).^2) <= radius^2;

        %Create a circle cutout
        cutout = image;
        for channel = 1:size(image, 3)
            cutout(:,:,channel) = image(:,:,channel) .* uint8(mask_circle);
        end
        
        
        %Get a zoomed cutout without black borders using bounding boxes
        boundingBox = regionprops(mask_circle, 'BoundingBox');
        boundingBox = boundingBox.BoundingBox;

        %Final cut-out of well
        cutout_final = imcrop(cutout, boundingBox);

        %Suplot of important steps of algorithm 
        figure(2);
        subplot(1,3,1)
        imshow(cutout_final);
        title('1');

        %Detecting the dimensions of the cropped image and the number of black pixels representing the edges
        background_of_cutout = sum(cutout_final(:) == 0); 
        [rows, cols, ~] = size(cutout_final); 
        number_of_pixels_cutout = rows * cols;

        %Dimensioning to detect colonies well_areating more than 50 cells 
        sum_rows = sum(mask_circle, 2);
        sum_columns = sum(mask_circle,1);
        pixels_diameter = max(sum_rows);
        pixels_diameter_2 = max(sum_columns);
        well_area = well_area;
        diameter_mm = (sqrt(well_area/pi)*2)*10;
        pixel = diameter_mm/pixels_diameter_2;
        number_of_pixels = colony_minimum/(pixel*pixel);
        number_of_pixels = ceil(number_of_pixels);

        %Kmeans segmentation 
        img_rgb = cutout_final;
        img_double = double(reshape(img_rgb, [], 3));

        %Splitting the image in half
        [rows, cols, ~] = size(img_rgb);
        half_rows = round(rows / 2);
        half_cols = round(cols / 2);

        %Number_of_clusters
        num_clusters = 3;

        %Inicialization of the result segmented image
        segmentedImg = zeros(rows, cols);

        %For loop for each half
        for j = 1:2
            %Actual half
            row_start = half_rows * (j - 1) + 1;
            row_end = min(half_rows * j, rows);

            %Actual half image
            img_half = img_rgb(: ,row_start:row_end, :);

            %Kmeans on half image
            img_half_double = double(reshape(img_half, [], 3));
            [~, centers] = kmeans(img_half_double, num_clusters);
            [~, cluster_idx] = min(pdist2(img_half_double, centers), [], 2);
            segmentedImg(:,row_start:row_end, :) = reshape(cluster_idx, size(img_half, 1), size(img_half, 2), []);

            %Number of pixels in each cluster
            num_pixels_in_clusters = histcounts(segmentedImg, 1:(num_clusters + 1));

            %Find the smallest cluster
            [~, min_cluster_idx] = min(num_pixels_in_clusters);

            %Create binary mask for smallest cluster
            min_cluster_binary = zeros(size(segmentedImg));
            min_cluster_binary(segmentedImg == min_cluster_idx) = 1;

        end


        %Suplot of important steps of algorithm 
        figure(2)
        subplot(1,3,2)
        imshow(min_cluster_binary);
        title('2');
        
        %Application of the Watershed algorithm to separate individual colonies
        binaryImage = ~bwareaopen(~min_cluster_binary, 20);
        distanceTransform = -bwdist(~binaryImage);
        watershedLabels = watershed(distanceTransform);

        %Since this is an oversegmentation, the following adjustment
        binaryImage2 = binaryImage;
        binaryImage2(watershedLabels == 0) = 0;

        maskMinima = imextendedmin(distanceTransform,2); 

        imposedMinima = imimposemin(distanceTransform,maskMinima); 
        watershedLabels2 = watershed(imposedMinima);
        resultImage = binaryImage;
        
        %Final image after Watershed algorithm
        resultImage(watershedLabels2 == 0) = 0; 

        
        %Call function for calculation colonies
        [valid_colony, count_valid_colony] = Count_Colonies(resultImage,number_of_pixels);

        %Colour separation of true and false colonies
        [rows, columns, ~] = size(valid_colony);

        new_image = zeros(rows, columns, 3, 'uint8'); 
        for i = 1:rows
            for j = 1:columns
                %Valid colonies will be red
                if valid_colony(i, j) > 0
                    new_image(i, j, :) = [255, 0, 0];
                end
            end
        end

        binaryImageRGB = im2uint8(cat(3, resultImage, resultImage, resultImage));
        combinedImage = imfuse(binaryImageRGB, new_image,'blend'); 

        %Extracted color channels
        redChannel = combinedImage(:, :, 1);
        greenChannel = combinedImage(:, :, 2);
        blueChannel = combinedImage(:, :, 3);

        %Find gray pixels
        blackPixels = redChannel == 128 & greenChannel  == 128 & blueChannel  == 128;
        redPixels = redChannel == 255 & greenChannel  == 128 & blueChannel  == 128;
        number_of_nonvalid_pixels = sum(blackPixels(:) == 1);
        number_of_valid_pixels = sum(redPixels(:) == 1);

        %Gray pixels coloring by cian (false colonies)
        redChannel(blackPixels) = 0;
        greenChannel(blackPixels) = 255;
        blueChannel(blackPixels) = 255;

        %Kombination of images and final one for color separation of
        %colonies
        rgbImage = cat(3, redChannel, greenChannel, blueChannel);
        subplot(1,3,3)
        imshow(rgbImage);
        title('3');

        %Calculation of colony area
        ratio = number_of_valid_pixels/(number_of_pixels_cutout - background_of_cutout);
        ratio_nonvalid = number_of_nonvalid_pixels/(number_of_pixels_cutout - background_of_cutout);
        ratio_well = (number_of_pixels_cutout - background_of_cutout-number_of_valid_pixels-number_of_nonvalid_pixels)/(number_of_pixels_cutout - background_of_cutout);

        disp(['Ratio of pixels of interest to colony in the well: ', num2str(ratio)]);
        disp(['Ratio of pixels of interest to nonvalid colony in the well: ', num2str(ratio_nonvalid)]);
        disp(['Ratio of pixels, which forms well: ', num2str(ratio_well)]);
        
        
        %Saving results
        area_valid_cm2 = ratio*well_area;
        area_nonvalid_cm2 = ratio_nonvalid*well_area;
        area_well_cm2 = ratio_well*well_area;
        disp(['Colony area: ', num2str(area_valid_cm2), ' cm2']);
        disp(['Area nonvalid colony: ', num2str(area_nonvalid_cm2), ' cm2']);
        disp(['Area well: ', num2str(area_well_cm2), ' cm2']);
        results(i+1,2)={area_valid_cm2};
        results(i+1,3)={count_valid_colony};
    end
    
end


%Function for finding valid colonies
function [valid_image, count_valid_colony]  = Count_Colonies(binaryImage, minimalSize)

    %Find connected area on image
    connected_component = bwconncomp(binaryImage);

    %Get colony property
    Colony_property = regionprops(connected_component, 'Centroid', 'Area', 'PixelList');

    %Filtration colony by minimal size
    valid_colony = Colony_property([Colony_property.Area] >= minimalSize);

    %Binary image with only walid colonies
    valid_image = false(size(binaryImage));

    for i = 1:length(valid_colony)
        pixelList = valid_colony(i).PixelList;
        indices = sub2ind(size(binaryImage), pixelList(:, 2), pixelList(:, 1));
        valid_image(indices) = true;
    end
       
    %Count of valid colonies
    count_valid_colony = length(valid_colony);
    %for i = 1:count_valid_colony
    %    text(valid_colony(i).Centroid(1), valid_colony(i).Centroid(2), num2str(i), 'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold');
    %end

    fprintf('Count of valid colonies (size >= %d): %d\n', minimalSize, count_valid_colony);

end



