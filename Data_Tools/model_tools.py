# imports


#  You can use this function to see how your model improves over time
def plot_history(history, metrics=None, ylim=None):
    plt.plot(history.history['loss'], label='training')
    plt.plot(history.history['val_loss'], label='testing')
    plt.title('Loss')
    if ylim:
        plt.ylim(ylim)
    plt.legend()
    plt.show()
    if metrics:
        for metric in metrics:
            plt.plot(history.history[metric], label=f'training {metric}')
            plt.plot(history.history[f'val_{metric}'], label=f'testing {metric}')
            plt.legend()
            plt.title(metric)
            if ylim:
                plt.ylim(ylim)
            plt.show()