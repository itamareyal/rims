'''
video_gen.py

generates a video from a collection of output histogram over time.
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import os
import cv2
from PIL import Image


'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''


def generate_video_from_frames(path_to_frames, title):
    # Description: Resize imgs to fit video.
    # Parameters: path_to_frames: location of all histograms to become the video
    #             title: name of the video created
    # Return: calls the video generation function.

    os.chdir(path_to_frames)
    mean_height = 0
    mean_width = 0
    num_of_images = len(os.listdir('.'))

    for file in os.listdir('.'):
        im = Image.open(file)
        width, height = im.size
        mean_width += width
        mean_height += height

    mean_width = int(mean_width / num_of_images)
    mean_height = int(mean_height / num_of_images)

    for file in os.listdir('.'):
        if file.endswith(".jpg") or file.endswith(".jpeg") or file.endswith("png") or file.endswith("JPEG"):
            # opening image using PIL Image
            im = Image.open(file)

            # resizing
            imResize = im.resize((mean_width, mean_height), Image.ANTIALIAS)
            imResize.save(file, 'JPEG', quality=95)

    # Calling the generate_video function
    generate_video(path_to_frames, title)


def generate_video(path_to_frames, title):
    image_folder = '.'  # make sure to use your folder
    video_name = title

    images = [img for img in os.listdir(image_folder)
              if img.endswith(".jpg") or
              img.endswith(".jpeg") or
              img.endswith(".JPEG") or
              img.endswith("png")]

    images = sorted(images, key=sort_by_title)
    # Array images should only consider
    # the image files ignoring others if any
    print("Video saved to "+path_to_frames + " as "+title)

    frame = cv2.imread(os.path.join(image_folder, images[0]))

    # setting the frame width, height width
    # the width, height of first image
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name, 0, 5, (width, height))

    # Appending the images to the video one by one
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

        # Deallocating memories taken for window creation
    cv2.destroyAllWindows()
    video.release()  # releasing the video generated


def sort_by_title(x):
    num = x.split('.')[0]
    num = num.split('_')[1]
    return int(num)

#generate_video_from_frames('RIMS output plots\\01-09-21_003124\\frames','clip2.avi')