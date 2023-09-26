import os
import boto3
import tarfile
import pandas as pd
import glob
import csv


def download_tar_files_from_directory(bucket_name, directory_name, destination_folder):
    s3 = boto3.client('s3')

    # List all objects in the directory
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=directory_name)

    # Iterate over the objects and download tar files
    for obj in response.get('Contents', []):
        key = obj['Key']
        filename = os.path.basename(key)

        # Check if the object is a tar file
        if filename.endswith('.tar'):
            # Download the tar file
            destination_path = os.path.join(destination_folder, filename)
            if not os.path.exists(destination_path):
                print(f'Downloading {filename}')
                s3.download_file(bucket_name, key, destination_path)
            

def unpack_tar_files(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # only unpack if it's a tar file and not already unpacked
        if filename.endswith('.tar'):
            print(f'Unpacking {filename}')
            with tarfile.open(file_path, 'r') as tar:
                tar.extractall(directory)
                print(f'Unpacked {filename}')
            os.remove(file_path)