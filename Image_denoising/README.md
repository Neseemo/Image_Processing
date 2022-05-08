{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "name": "README.md",
      "provenance": []
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "source": [
        "# Image Denoising\n",
        "In this notebook are presented two algorithms for image denoising:\n",
        "\n",
        "\n",
        "*   Hard thresholding of DCT representation of the image patches \n",
        "*   smoothing by convolution (average of the pixels)\n",
        "\n",
        "The first one works better but assumes the variance of the added noise to be known (wich is not obviously the case in reality). However, in the final part of the notebook a noise estimation technique is presented and, as shown, the estimated variance is very close to the actual one so the first algorith remains the best.\n"
      ],
      "metadata": {
        "id": "s5GKu41W708e"
      }
    }
  ]
}