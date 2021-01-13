import { DomainItem } from './Item'
import { v4 as uuidv4 } from 'uuid'

export const DomainList = props => (
  <div>
    <div>
      <button type="button" onClick={props.handleAdd}>Add</button>
    </div>
    <ul>
      {props.domains.map((domain, index) => {
        if (!domain.uuid)
          domain.uuid = uuidv4()
        return (
          <DomainItem
            key={domain.uuid}
            data={domain}
            handleRemove={props.handleRemove(index)}
            handleChange={props.handleChange(index)}
          />
        )
      })}
    </ul>
  </div>
)
