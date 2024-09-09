#pragma once
#include <iostream>
#include <librdkafka/rdkafka.h>
#include <string>

class KafkaProducer
{
  public:
    KafkaProducer(const std::string& brokers, const std::string& topic)
    {
        rd_kafka_conf_t* conf = rd_kafka_conf_new();

        char errstr[512];
        if (rd_kafka_conf_set(conf, "bootstrap.servers", brokers.c_str(), errstr, sizeof(errstr)) != RD_KAFKA_CONF_OK)
        {
            std::cerr << "Error setting brokers: " << errstr << std::endl;
            throw std::runtime_error("Failed to configure Kafka");
        }

        producer = rd_kafka_new(RD_KAFKA_PRODUCER, conf, errstr, sizeof(errstr));
        if (!producer)
        {
            std::cerr << "Failed to create producer: " << errstr << std::endl;
            throw std::runtime_error("Failed to create Kafka producer");
        }

        this->topic = rd_kafka_topic_new(producer, topic.c_str(), nullptr);
        if (!this->topic)
        {
            std::cerr << "Failed to create topic: " << errstr << std::endl;
            throw std::runtime_error("Failed to create Kafka topic");
        }
    }

    ~KafkaProducer()
    {
        rd_kafka_flush(producer, 10 * 1000);
        rd_kafka_topic_destroy(topic);
        rd_kafka_destroy(producer);
    }

    void send(const std::string& message)
    {
        rd_kafka_produce(
            topic,
            RD_KAFKA_PARTITION_UA,
            RD_KAFKA_MSG_F_COPY,
            const_cast<char*>(message.c_str()),
            message.size(),
            nullptr,
            0,
            nullptr);
    }

  private:
    rd_kafka_t* producer;
    rd_kafka_topic_t* topic;
};
