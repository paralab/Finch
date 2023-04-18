//
// Created by milinda on 11/19/18.
//

/**
 * @brief Simple class to manage async data transfer in the ODA class.
 * */


#ifndef DENDRO_KT_UPDATECTX_H
#define DENDRO_KT_UPDATECTX_H

#include "nsort.h"

#include "mpi.h"
/// #include <vector>
#include <memory>
#include <typeinfo>

namespace ot {

    class AsyncExchangeContex {

        private :
            size_t m_bufferType;

            /** pointer to the variable which perform the ghost exchange */
            void* m_uiBuffer;

            /** pointer to the send buffer*/
            void* m_uiSendBuf;

            /** pointer to the send buffer*/
            void* m_uiRecvBuf;

            std::unique_ptr<MPI_Request[]> m_uiUpstRequests;
            std::unique_ptr<MPI_Request[]> m_uiDnstRequests;

        public:
            /**@brief creates an async ghost exchange contex*/
            AsyncExchangeContex(const void* var, size_t bufferType, unsigned int nUpstProcs, unsigned int nDnstProcs)
              : m_uiUpstRequests{new MPI_Request[nUpstProcs]},
                m_uiDnstRequests{new MPI_Request[nDnstProcs]}
            {
                m_bufferType = bufferType;
                m_uiBuffer=(void*)var;
                m_uiSendBuf=NULL;
                m_uiRecvBuf=NULL;
            }

            /**@brief Move constructor. */
            AsyncExchangeContex(AsyncExchangeContex &&that)
              :  m_uiUpstRequests{std::move(that.m_uiUpstRequests)},
                 m_uiDnstRequests{std::move(that.m_uiDnstRequests)},
                 m_bufferType{that.m_bufferType},
                 m_uiBuffer{that.m_uiBuffer},
                 m_uiSendBuf{that.m_uiSendBuf},
                 m_uiRecvBuf{that.m_uiRecvBuf}
            {}

            /**@brief Move assignment operator. */
            void operator= (AsyncExchangeContex &&that)
            {
              m_uiUpstRequests = std::move(that.m_uiUpstRequests);
              m_uiDnstRequests = std::move(that.m_uiDnstRequests);
              m_bufferType = that.m_bufferType;
              m_uiBuffer = that.m_uiBuffer;
              m_uiSendBuf = that.m_uiSendBuf;
              m_uiRecvBuf = that.m_uiRecvBuf;
            }

            /**@brief allocates send buffer for ghost exchange*/
            inline void allocateSendBuffer(size_t bytes)
            {
                m_uiSendBuf=malloc(bytes);
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void allocateRecvBuffer(size_t bytes)
            {
                m_uiRecvBuf=malloc(bytes);
            }

            /**@brief allocates send buffer for ghost exchange*/
            inline void deAllocateSendBuffer()
            {
                free(m_uiSendBuf);
                m_uiSendBuf=NULL;
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void deAllocateRecvBuffer()
            {
                free(m_uiRecvBuf);
                m_uiRecvBuf=NULL;
            }

            inline size_t getBufferType() { return m_bufferType; }

            inline void* getSendBuffer() { return m_uiSendBuf;}
            inline void* getRecvBuffer() { return m_uiRecvBuf;}

            inline void* getBuffer() const {return m_uiBuffer;}

            /** @note Remember not to copy-assign MPI_Request. */
            inline MPI_Request * getUpstRequestList(){ return m_uiUpstRequests.get();}
            inline MPI_Request * getDnstRequestList(){ return m_uiDnstRequests.get();}

            bool operator== (const AsyncExchangeContex &other) const{
                return( m_uiBuffer == other.m_uiBuffer );
            }

            ~AsyncExchangeContex() {
            }

    };

} //end namespace

#endif //DENDRO_KT_UPDATECTX_H
